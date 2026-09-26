#include "mpc_cost.h"

#include <stddef.h>

/* ============================================================
 * mpc_cost.c
 *
 * MPC cost function reformulated as QP.
 *
 * MATLAB source tracking:
 *   kalmanCoastingPhase.m
 * ============================================================ */

/* ============================================================
 * Private helper functions
 * ============================================================ */

/* ------------------------------------------------------------
 * Build x vector from state estimate
 *
 * MATLAB tracking:
 *   x = [z; v]
 * ------------------------------------------------------------ */
static void build_state_vector(const MPC_StateEstimate *y_estim, double x[MPC_N_STATES]);

/* ------------------------------------------------------------
 * Build diagonal Q vector per stacked output component
 *
 * MATLAB tracking:
 *   Q_step = [q_z, q_v];
 *   Q = diag(repmat(Q_step, 1, Np));
 *
 * Since Q is diagonal, we store only its diagonal.
 * ------------------------------------------------------------ */
static void build_Q_diag(const MPC_Weights *weights, double Q_diag[MPC_RP_SIZE]);

/* ------------------------------------------------------------
 * Build diagonal R vector
 *
 * MATLAB tracking:
 *   R = 1
 *
 * More generally for delta_U stacked over horizon:
 *   R = r_du * I
 * ------------------------------------------------------------ */
static void build_R_diag(const MPC_Weights *weights, double R_diag[MPC_DELTA_U_SIZE]);

/* ------------------------------------------------------------
 * Compute affine free response:
 *
 *   Y0 = Phi_extend*x + Phiref*REF
 *
 * MATLAB tracking:
 *   Y = Phi_extend*x + Gamma_extend*delta_U + Phiref*REF
 * ------------------------------------------------------------ */
static void build_free_response(
    const double Phi_extend[MPC_RP_SIZE][MPC_N_STATES],
    const double x[MPC_N_STATES],
    const double Phiref[MPC_RP_SIZE][MPC_N_STATES],
    const double REF[MPC_N_STATES],
    double Y0[MPC_RP_SIZE]
);

/* ============================================================
 * Public functions
 * ============================================================ */

int MPC_objfun_to_qp(const MPC_StateEstimate *y_estim, const double Phi_extend[MPC_RP_SIZE][MPC_N_STATES], const double Gamma_extend[MPC_RP_SIZE][MPC_DELTA_U_SIZE], const double r_p[MPC_RP_SIZE], const double Phiref[MPC_RP_SIZE][MPC_N_STATES], const double REF[MPC_N_STATES], const MPC_Weights *weights, double P[MPC_DELTA_U_SIZE][MPC_DELTA_U_SIZE], double q[MPC_DELTA_U_SIZE])
{
    double x[MPC_N_STATES];
    double Q_diag[MPC_RP_SIZE];
    double R_diag[MPC_DELTA_U_SIZE];
    double Y0[MPC_RP_SIZE];
    double e0[MPC_RP_SIZE];
    int i, j, k;

    if (y_estim == NULL || Phi_extend == NULL || Gamma_extend == NULL ||
        r_p == NULL || Phiref == NULL || REF == NULL ||
        weights == NULL || P == NULL || q == NULL) {
        return -1;
    }

    build_state_vector(y_estim, x);
    build_Q_diag(weights, Q_diag);
    build_R_diag(weights, R_diag);
    build_free_response(Phi_extend, x, Phiref, REF, Y0);

    /* ========================================================
     * MATLAB block:
     *   Y = Phi_extend*x + Gamma_extend*delta_U + Phiref*REF;
     *
     * Define:
     *   Y = Y0 + Gamma_extend * delta_U
     *
     * where:
     *   Y0 = Phi_extend*x + Phiref*REF
     * ======================================================== */

    /* ========================================================
     * MATLAB block:
     *   MPC_objfun = 1/2*(r_p-Y)'*Q*(r_p-Y) + 1/2*delta_U'*R*delta_U;
     *
     * Let:
     *   e0 = Y0 - r_p
     *
     * Then:
     *   J(delta_U) = 1/2 * (Gamma*du + e0)' Q (Gamma*du + e0)
     *                + 1/2 * du' R du
     *
     * Expanding:
     *   J(delta_U) = 1/2 * du' * (Gamma'QGamma + R) * du
     *                + (Gamma'Qe0)' * du
     *                + constant
     *
     * Therefore for QP form:
     *   1/2 * du' P du + q' du
     *
     * we have:
     *   P = Gamma'QGamma + R
     *   q = Gamma'Qe0
     * ======================================================== */

    for (i = 0; i < MPC_RP_SIZE; ++i) {
        e0[i] = Y0[i] - r_p[i];
    }

    /* Build P = Gamma'QGamma + R */
    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        for (j = 0; j < MPC_DELTA_U_SIZE; ++j) {
            double sum = 0.0;

            for (k = 0; k < MPC_RP_SIZE; ++k) {
                sum += Gamma_extend[k][i] * Q_diag[k] * Gamma_extend[k][j];
            }

            if (i == j) {
                sum += R_diag[i];
            }

            P[i][j] = sum;
        }
    }

    /* Build q = Gamma'Qe0 */
    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        double sum = 0.0;

        for (k = 0; k < MPC_RP_SIZE; ++k) {
            sum += Gamma_extend[k][i] * Q_diag[k] * e0[k];
        }

        q[i] = sum;
    }

    return 0;
}

/* ============================================================
 * Private helper implementations
 * ============================================================ */

static void build_state_vector(const MPC_StateEstimate *y_estim, double x[MPC_N_STATES])
{
    x[0] = y_estim->z;
    x[1] = y_estim->v;
}

static void build_Q_diag(const MPC_Weights *weights, double Q_diag[MPC_RP_SIZE])
{
    int i;

    /* MATLAB tracking:
     *   Q_step = [q_z, q_v];
     *   Q = diag(repmat(Q_step, 1, Np));
     */
    for (i = 0; i < MPC_NP; ++i) {
        Q_diag[2 * i + 0] = weights->q_z;
        Q_diag[2 * i + 1] = weights->q_v;
    }
}

static void build_R_diag(const MPC_Weights *weights, double R_diag[MPC_DELTA_U_SIZE])
{
    int i;

    /* MATLAB tracking:
     *   R = r_du * I
     */
    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        R_diag[i] = weights->r_du;
    }
}

static void build_free_response(const double Phi_extend[MPC_RP_SIZE][MPC_N_STATES], const double x[MPC_N_STATES], const double Phiref[MPC_RP_SIZE][MPC_N_STATES], const double REF[MPC_N_STATES], double Y0[MPC_RP_SIZE])
{
    int i;

    for (i = 0; i < MPC_RP_SIZE; ++i) {
        Y0[i] =
            Phi_extend[i][0] * x[0] +
            Phi_extend[i][1] * x[1] +
            Phiref[i][0] * REF[0] +
            Phiref[i][1] * REF[1];
    }
}