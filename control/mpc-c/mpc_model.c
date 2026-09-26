#include "mpc_model.h"

#include <math.h>
#include <stddef.h>

/* ============================================================
 * mpc_model.c
 *
 * Model-related functions used by the controller.
 *
 * MATLAB source tracking:
 *   kalmanCoastingPhase.m
 *
 * This file groups the model / aerodynamic / local linearization
 * functions needed by the MPC block.
 * ============================================================ */

/* ============================================================
 * Private helper functions
 * ============================================================ */

/* ------------------------------------------------------------
 * Linear interpolation helper
 *
 * MATLAB tracking:
 *   interp1(...)
 *
 * Assumed table format:
 *   table[i][0] = y value
 *   table[i][1] = x value
 *
 * This matches:
 *   interp1(air_dens_alt(:,2), air_dens_alt(:,1), altitude)
 *   interp1(CD_speed(:,2), CD_speed(:,1), velocity)
 * ------------------------------------------------------------ */
static int interp1_table_col01(const double (*table)[2], int len, double xq, double *yq);

/* ------------------------------------------------------------
 * Matrix exponential for a 4x4 matrix using truncated series.
 *
 * Used for gamma_integral exact MATLAB-style implementation.
 *
 * MATLAB tracking:
 *   EM = expm(M * h);
 *
 * Note:
 * This is a compact implementation for small fixed size only.
 * ------------------------------------------------------------ */
static void mat4_identity(double A[4][4]);
static void mat4_mul(const double A[4][4], const double B[4][4], double C[4][4]);
static void mat4_add_inplace(double A[4][4], const double B[4][4]);
static void mat4_scale(const double A[4][4], double s, double B[4][4]);
static void mat4_expm_series(const double A[4][4], double E[4][4]);

/* ============================================================
 * Public functions
 * ============================================================ */

int compute_current_parameters(const double (*air_dens_alt)[2], int air_dens_alt_len, const double (*CD_speed)[2], int cd_speed_len, double altitude, double velocity, double *air_dens, double *drag_coeff)
{
    int status;

    if (air_dens_alt == NULL || CD_speed == NULL || air_dens == NULL || drag_coeff == NULL) {
        return -1;
    }

    /* ========================================================
     * MATLAB block:
     *   air_dens = interp1(air_dens_alt(:,2), air_dens_alt(:,1), y_estim(1,s));
     * ======================================================== */
    status = interp1_table_col01(air_dens_alt, air_dens_alt_len, altitude, air_dens);
    if (status != 0) {
        return -2;
    }

    /* ========================================================
     * MATLAB block:
     *   drag_coeff = interp1(CD_speed(:,2), CD_speed(:,1), y_estim(2,s));
     *
     * IMPORTANT:
     * This follows your MATLAB literally: signed velocity is used,
     * not abs(velocity).
     * ======================================================== */
    status = interp1_table_col01(CD_speed, cd_speed_len, velocity, drag_coeff);
    if (status != 0) {
        return -3;
    }

    return 0;
}

int lin_system(double Cd_global, double rho_der, double S, double mass, double air_dens, const MPC_StateEstimate *y_estim, double Sab, double ab_drag_coeff, double A_c[2][2], double B_c[2])
{
    double v;

    if (y_estim == NULL || A_c == NULL || B_c == NULL) {
        return -1;
    }

    v = y_estim->v;

    /* ========================================================
     * MATLAB block:
     *
     * function [A_c, B_c] = lin_system(Cd_global, rho_der, S, mass,
     *                                  air_dens, x, Sab, Cd_ab)
     *
     * A_c = [0, 1;
     *   -0.5/mass*rho_der*x(2)^2*S*Cd_global,
     *   -1/mass*air_dens*x(2)*S*Cd_global];
     *
     * B_c = [0;
     *   -0.5/mass*air_dens*x(2)^2*Sab*Cd_ab];
     * ======================================================== */

    

    A_c[0][0] = 0.0;
    A_c[0][1] = 1.0;
    A_c[1][0] = -0.5 / mass * rho_der * v * v * S * Cd_global;
    A_c[1][1] = -1.0 / mass * air_dens * v * S * Cd_global;

    B_c[0] = 0.0;
    B_c[1] = -0.5 / mass * air_dens * v * v * Sab * ab_drag_coeff;

    return 0;
}

int gamma_integral(const double A[2][2], double h, double Gamma_hat[2][2])
{
    double M[4][4] = {{0.0}};
    double Mh[4][4] = {{0.0}};
    double EM[4][4] = {{0.0}};
    int i, j;

    if (A == NULL || Gamma_hat == NULL) {
        return -1;
    }

    /* ========================================================
     * MATLAB block:
     *
     * function integral_gamma = gamma_integral(A, h)
     *     n = size(A,1);
     *     M = [A, eye(n);
     *          zeros(n), zeros(n)];
     *     EM = expm(M * h);
     *     integral_gamma = EM(1:n, n+1:end);
     * end
     * ======================================================== */

    /* M = [A, I;
     *      0, 0]
     */
    M[0][0] = A[0][0];
    M[0][1] = A[0][1];
    M[0][2] = 1.0;
    M[0][3] = 0.0;

    M[1][0] = A[1][0];
    M[1][1] = A[1][1];
    M[1][2] = 0.0;
    M[1][3] = 1.0;

    M[2][0] = 0.0;
    M[2][1] = 0.0;
    M[2][2] = 0.0;
    M[2][3] = 0.0;

    M[3][0] = 0.0;
    M[3][1] = 0.0;
    M[3][2] = 0.0;
    M[3][3] = 0.0;

    /* Mh = M * h */
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            Mh[i][j] = M[i][j] * h;
        }
    }

    /* EM = expm(M*h) */
    mat4_expm_series(Mh, EM);

    /* integral_gamma = EM(1:n, n+1:end) */
    Gamma_hat[0][0] = EM[0][2];
    Gamma_hat[0][1] = EM[0][3];
    Gamma_hat[1][0] = EM[1][2];
    Gamma_hat[1][1] = EM[1][3];

    return 0;
}

int compute_dynamics(const MPC_StateEstimate *y_estim, double U, double air_dens, double drag_coeff, double ab_drag_coeff, double Fref[2])
{
    double speed;
    double deflection;
    double acc;

    if (y_estim == NULL || Fref == NULL) {
        return -1;
    }

    speed = y_estim->v;
    deflection = U;

    /* ========================================================
     * MATLAB block:
     *
     * function [dXdt] = compute_dynamics(X, U, air_dens, CD, CDab)
     *     altitude = X(1);
     *     speed = X(2);
     *     deflection = U;
     *
     *     acc = -9.77 - 0.5 * air_dens * 9.5e-3 * speed^2 / 16.4813 ...
     *           * (CD + CDab * 4e-3/9.5e-3 * deflection);
     *
     *     dXdt = [speed; acc];
     * end
     * ======================================================== */

    acc = -9.77 - 0.5 * air_dens * 9.5e-3 * speed * speed / 16.4813
          * (drag_coeff + ab_drag_coeff * 4e-3 / 9.5e-3 * deflection);

    Fref[0] = speed;
    Fref[1] = acc;

    return 0;
}

/* ============================================================
 * Private helper implementations
 * ============================================================ */

static int interp1_table_col01(const double (*table)[2], int len, double xq, double *yq)
{
    int i;

    if (table == NULL || yq == NULL || len < 2) {
        return -1;
    }

    /* Saturation on lower bound */
    if (xq <= table[0][1]) {
        *yq = table[0][0];
        return 0;
    }

    /* Saturation on upper bound */
    if (xq >= table[len - 1][1]) {
        *yq = table[len - 1][0];
        return 0;
    }

    /* Search interval */
    for (i = 0; i < len - 1; ++i) {
        double x0 = table[i][1];
        double x1 = table[i + 1][1];

        if (xq >= x0 && xq <= x1) {
            double y0 = table[i][0];
            double y1 = table[i + 1][0];
            double alpha;

            if (x1 == x0) {
                *yq = y0;
                return 0;
            }

            alpha = (xq - x0) / (x1 - x0);
            *yq = y0 + alpha * (y1 - y0);
            return 0;
        }
    }

    return -2;
}

static void mat4_identity(double A[4][4])
{
    int i, j;

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            A[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

static void mat4_mul(const double A[4][4], const double B[4][4], double C[4][4])
{
    int i, j, k;

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            double sum = 0.0;
            for (k = 0; k < 4; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

static void mat4_add_inplace(double A[4][4], const double B[4][4])
{
    int i, j;

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            A[i][j] += B[i][j];
        }
    }
}

static void mat4_scale(const double A[4][4], double s, double B[4][4])
{
    int i, j;

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            B[i][j] = s * A[i][j];
        }
    }
}

static void mat4_expm_series(const double A[4][4], double E[4][4])
{
    double term[4][4];
    double temp[4][4];
    int k;
    double inv_k = 1.0;

    /* E = I */
    mat4_identity(E);

    /* term = I, then recurrence:
     * term_k = term_{k-1} * A / k
     */
    mat4_identity(term);

    /* Truncated series for small fixed-size matrix.
     * Good enough as a first implementation for moderate h.
     */
    for (k = 1; k <= 20; ++k) {
        mat4_mul(term, A, temp);
        inv_k = 1.0 / (double)k;
        mat4_scale(temp, inv_k, term);
        mat4_add_inplace(E, term);
    }
}