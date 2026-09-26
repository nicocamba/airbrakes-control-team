#include "mpc_controller.h"

#include <stddef.h>
#include <string.h>

#include "mpc_model.h"
#include "mpc_prediction.h"
#include "mpc_constraints.h"
#include "mpc_cost.h"
#include "mpc_solver.h"

/* ============================================================
 * mpc_controller.c
 *
 * High-level MPC controller sequence.
 *
 * MATLAB source tracking:
 *   kalmanCoastingPhase.m
 *
 * This file keeps the top-level control flow:
 *   reference -> parameters -> linearization -> discretization
 *   -> prediction -> constraints -> cost -> solver -> control u
 * ============================================================ */

/* ============================================================
 * Private helper functions (local to this .c only)
 * ============================================================ */

/* ------------------------------------------------------------
 * Build reference horizon r_p
 *
 * MATLAB tracking:
 *   for i = 1:Np
 *       r_p(n*i-1:n*i) = [altitude_ref(s+i-1); velocity_ref(s+i-1)];
 *   end
 *
 * C convention:
 * - s is 0-based
 * - i goes 0..Np-1
 * - r_p[2*i + 0] = altitude_ref[s+i]
 * - r_p[2*i + 1] = velocity_ref[s+i]
 * ------------------------------------------------------------ */
static int build_reference_horizon(const double *altitude_ref, const double *velocity_ref, int ref_length, int s, double r_p[MPC_RP_SIZE]);

/* ------------------------------------------------------------
 * Warm start initialization for delta_U0
 *
 * MATLAB tracking:
 *   delta_U0 = zeros(m*Np,1);
 *   for i=1:Np-1
 *       delta_U0(i) = delta_Uprev(i+1);
 *   end
 *
 * Since m = 1:
 * - shift previous optimal sequence to the left
 * - set the last element to zero
 * ------------------------------------------------------------ */
static void build_warm_start(const double delta_Uprev[MPC_DELTA_U_SIZE], double delta_U0[MPC_DELTA_U_SIZE]);

/* ------------------------------------------------------------
 * First-order approximation of Phi
 *
 * MATLAB tracking:
 *   Phi = expm(A_c*h);
 *
 * Current scaffold:
 *   Phi ≈ I + A_c*h
 *
 * Note:
 * This is only a temporary approximation for structure/testing.
 * ------------------------------------------------------------ */
static void compute_Phi_first_order(const double A_c[2][2], double h, double Phi[2][2]);

/* ============================================================
 * Public functions
 * ============================================================ */

void mpc_init_delta_u(double delta_U[MPC_DELTA_U_SIZE])
{
    int i;

    if (delta_U == NULL) {
        return;
    }

    /* MATLAB equivalent:
     *   delta_U = zeros(m*Np,1);
     */
    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        delta_U[i] = 0.0;
    }
}

int mpc_controller_step(const MPC_StateEstimate *y_estim, const double *altitude_ref, const double *velocity_ref, int ref_length, int s, const double (*air_dens_alt)[2], int air_dens_alt_len, const double (*CD_speed)[2], int cd_speed_len, double uprev, const double delta_Uprev[MPC_DELTA_U_SIZE], const MPC_Params *params, const MPC_Weights *weights, const MPC_Bounds *bounds, MPC_Output *out)
{
    /* ========================================================
     * Local variables mirroring MATLAB names as much as possible
     * ======================================================== */

    double r_p[MPC_RP_SIZE];

    double air_dens = 0.0;
    double drag_coeff = 0.0;
    double Cd_global = 0.0;

    double A_c[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double B_c[2] = {0.0, 0.0};

    double Phi[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double Gamma_hat[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double Gamma[2] = {0.0, 0.0};

    double Fref[2] = {0.0, 0.0};

    double Phi_extend[MPC_RP_SIZE][MPC_N_STATES];
    double Gamma_extend[MPC_RP_SIZE][MPC_DELTA_U_SIZE];
    double Phiref[MPC_RP_SIZE][MPC_N_STATES];
    double REF[MPC_N_STATES];

    double A_cons[2 * MPC_NP][MPC_DELTA_U_SIZE];
    double b_cons[2 * MPC_NP];

    double P[MPC_DELTA_U_SIZE][MPC_DELTA_U_SIZE];
    double q[MPC_DELTA_U_SIZE];

    double delta_U0[MPC_DELTA_U_SIZE];
    double delta_U_opt[MPC_DELTA_U_SIZE];

    int i;
    int status;

    /* ========================================================
     * Input validation
     * ======================================================== */
    if (y_estim == NULL || altitude_ref == NULL || velocity_ref == NULL ||
        air_dens_alt == NULL || CD_speed == NULL ||
        params == NULL || weights == NULL || bounds == NULL || out == NULL) {
        return -1;
    }

    /* Safe output initialization */
    out->u = uprev;
    out->solver_status = -999;
    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        out->delta_U[i] = 0.0;
    }

    /* ========================================================
     * MATLAB block:
     *   for i = 1:Np
     *       r_p(n*i-1:n*i) = [altitude_ref(s+i-1); velocity_ref(s+i-1)];
     *   end
     * ======================================================== */
    status = build_reference_horizon(altitude_ref, velocity_ref, ref_length, s, r_p);
    if (status != 0) {
        return -2;
    }

    /* ========================================================
     * MATLAB block:
     *   [air_dens, drag_coeff] = compute_current_parameters(
     *       air_dens_alt, CD_speed, y_estim(1,s), y_estim(2,s));
     * ======================================================== */
    status = compute_current_parameters(air_dens_alt, air_dens_alt_len, CD_speed, cd_speed_len, y_estim->z, y_estim->v, &air_dens, &drag_coeff);
    if (status != 0) {
        return -3;
    }

    /* ========================================================
     * MATLAB block:
     *   Cd_global = drag_coeff + Sab / S * ab_drag_coeff*u;
     *
     * Here u is approximated with uprev for the linearization point.
     * ======================================================== */
    Cd_global = drag_coeff + (params->Sab / params->S) * params->ab_drag_coeff * uprev;

    /* ========================================================
     * MATLAB block:
     *   [A_c, B_c] = lin_system(Cd_global, rho_der, S, mass,
     *                           air_dens, y_estim(:,s),
     *                           Sab, ab_drag_coeff);
     * ======================================================== */
    status = lin_system(Cd_global, params->rho_der, params->S, params->mass, air_dens, y_estim, params->Sab, params->ab_drag_coeff, A_c, B_c);
    if (status != 0) {
        return -4;
    }

    /* ========================================================
     * MATLAB block:
     *   Phi = expm(A_c*h);
     *
     * Current scaffold:
     *   first-order approximation only
     * ======================================================== */
    compute_Phi_first_order(A_c, params->h, Phi);

    /* ========================================================
     * MATLAB block:
     *   Gamma_hat = gamma_integral(A_c, h);
     * ======================================================== */
    status = gamma_integral(A_c, params->h, Gamma_hat);
    if (status != 0) {
        return -5;
    }

    /* ========================================================
     * MATLAB block:
     *   Gamma = Gamma_hat * B_c;
     * ======================================================== */
    Gamma[0] = Gamma_hat[0][0] * B_c[0] + Gamma_hat[0][1] * B_c[1];
    Gamma[1] = Gamma_hat[1][0] * B_c[0] + Gamma_hat[1][1] * B_c[1];

    /* ========================================================
     * MATLAB block:
     *   Fref = compute_dynamics(y_estim(:,s), u, air_dens,
     *                           drag_coeff, ab_drag_coeff);
     *
     * Here u is approximated with uprev.
     * ======================================================== */
    status = compute_dynamics(y_estim, uprev, air_dens, drag_coeff, params->ab_drag_coeff, Fref);
    if (status != 0) {
        return -6;
    }

    /* ========================================================
     * MATLAB block:
     *   [Phi_extend, Gamma_extend, Phiref, REF] =
     *       prediction_matrix_MPC(n, m, Np, A_c, Phi, x,
     *                             Gamma_hat, Gamma, Fref);
     *
     * Current C architecture:
     * - prediction built around estimated state y_estim
     * ======================================================== */
    status = prediction_matrix_MPC(y_estim, A_c, Phi, Gamma_hat, Gamma, Fref, Phi_extend, Gamma_extend, Phiref, REF);
    if (status != 0) {
        return -7;
    }

    /* ========================================================
     * MATLAB block:
     *   delta_U0 = zeros(m*Np,1);
     *   for i=1:Np-1
     *       delta_U0(i) = delta_Uprev(i+1);
     *   end
     * ======================================================== */
    build_warm_start(delta_Uprev, delta_U0);

    /* ========================================================
     * MATLAB block:
     *   [A, b] = MPC_Lconsfun(Umin, Umax, uprev, E, H);
     *
     * Current C scaffold:
     * - scalar bounds are expanded internally across horizon
     * ======================================================== */
    status = MPC_Lconsfun(bounds->U_min, bounds->U_max, uprev, A_cons, b_cons);
    if (status != 0) {
        return -8;
    }

    /* ========================================================
     * MATLAB block:
     *   fun = @(delta_U) MPC_objfun(y_estim(1:2,s), Phi_extend,
     *       Gamma_extend, r_p, Q, R, delta_U, Phiref, REF);
     *
     * Current C/QP approach:
     * - build quadratic matrices P and q
     * ======================================================== */
    status = MPC_objfun_to_qp(y_estim, Phi_extend, Gamma_extend, r_p, Phiref, REF, weights, P, q);
    if (status != 0) {
        return -9;
    }

    /* ========================================================
     * MATLAB block:
     *   [delta_U,...] = fmincon(fun, delta_U0, A, b, [], [], [], [], [], options);
     *
     * Current C approach:
     * - solve equivalent QP through solver abstraction
     * ======================================================== */
    status = solve_qp(P, q, A_cons, b_cons, delta_U0, delta_U_opt);

    out->solver_status = status;

    if (status < 0) {
        out->u = uprev;
        for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
            out->delta_U[i] = 0.0;
        }
        return -10;
    }

    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        out->delta_U[i] = delta_U_opt[i];
    }

    /* ========================================================
     * MATLAB block:
     *   u = uprev + delta_U(1);
     *
     * C indexing:
     *   delta_U_opt[0]
     * ======================================================== */
    out->u = uprev + delta_U_opt[0];

    return 0;
}

/* ============================================================
 * Private helper implementations
 * ============================================================ */

static int build_reference_horizon(const double *altitude_ref, const double *velocity_ref, int ref_length, int s, double r_p[MPC_RP_SIZE])
{
    int i;

    if (altitude_ref == NULL || velocity_ref == NULL || r_p == NULL) {
        return -1;
    }

    for (i = 0; i < MPC_NP; ++i) {
        int idx_ref = s + i;
        int idx_rp = 2 * i;

        if (idx_ref < 0 || idx_ref >= ref_length) {
            return -2;
        }

        /* MATLAB equivalent:
         *   r_p(n*i-1:n*i) = [altitude_ref(s+i-1); velocity_ref(s+i-1)];
         */
        r_p[idx_rp + 0] = altitude_ref[idx_ref];
        r_p[idx_rp + 1] = velocity_ref[idx_ref];
    }

    return 0;
}

static void build_warm_start(const double delta_Uprev[MPC_DELTA_U_SIZE], double delta_U0[MPC_DELTA_U_SIZE])
{
    int i;

    if (delta_U0 == NULL) {
        return;
    }

    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        delta_U0[i] = 0.0;
    }

    if (delta_Uprev == NULL) {
        return;
    }

    for (i = 0; i < MPC_DELTA_U_SIZE - 1; ++i) {
        /* MATLAB equivalent:
         *   delta_U0(i) = delta_Uprev(i+1)
         * adapting MATLAB 1-based indexing to C 0-based indexing
         */
        delta_U0[i] = delta_Uprev[i + 1];
    }
}

static void compute_Phi_first_order(const double A_c[2][2], double h, double Phi[2][2])
{
    /* MATLAB target:
     *   Phi = expm(A_c*h);
     *
     * Current approximation:
     *   Phi ≈ I + A_c*h
     */

    Phi[0][0] = 1.0 + A_c[0][0] * h;
    Phi[0][1] =       A_c[0][1] * h;
    Phi[1][0] =       A_c[1][0] * h;
    Phi[1][1] = 1.0 + A_c[1][1] * h;
}