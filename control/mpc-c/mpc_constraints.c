#include "mpc_constraints.h"

#include <stddef.h>

/* ============================================================
 * mpc_constraints.c
 *
 * Linear inequality constraints for MPC.
 *
 * MATLAB source tracking:
 *   kalmanCoastingPhase.m
 * ============================================================ */

/* ============================================================
 * Private helper functions
 * ============================================================ */

/* ------------------------------------------------------------
 * Build E matrix
 *
 * MATLAB role:
 *   E * u
 *
 * For m = 1, E is simply a column vector of ones of size Np x 1.
 * ------------------------------------------------------------ */
static void build_E(double E[MPC_NP]);

/* ------------------------------------------------------------
 * Build H matrix
 *
 * MATLAB role:
 *   U = E*u + H*delta_U
 *
 * For m = 1, H is lower triangular with ones on and below diagonal.
 * ------------------------------------------------------------ */
static void build_H(double H[MPC_NP][MPC_DELTA_U_SIZE]);

/* ============================================================
 * Public functions
 * ============================================================ */

int MPC_Lconsfun(double U_min, double U_max, double uprev, double A_cons[2 * MPC_NP][MPC_DELTA_U_SIZE], double b_cons[2 * MPC_NP])
{
    double E[MPC_NP];
    double H[MPC_NP][MPC_DELTA_U_SIZE];
    int i, j;

    if (A_cons == NULL || b_cons == NULL) {
        return -1;
    }

    build_E(E);
    build_H(H);

    /* ========================================================
     * MATLAB block:
     *
     *   A = [-H;
     *         H];
     * ======================================================== */
    for (i = 0; i < MPC_NP; ++i) {
        for (j = 0; j < MPC_DELTA_U_SIZE; ++j) {
            A_cons[i][j] = -H[i][j];
            A_cons[MPC_NP + i][j] = H[i][j];
        }
    }

    /* ========================================================
     * MATLAB block:
     *
     *   b = [(-Umin + E*u);
     *         ( Umax - E*u)];
     * ======================================================== */
    for (i = 0; i < MPC_NP; ++i) {
        b_cons[i] = -U_min + E[i] * uprev;
        b_cons[MPC_NP + i] = U_max - E[i] * uprev;
    }

    return 0;
}

/* ============================================================
 * Private helper implementations
 * ============================================================ */

static void build_E(double E[MPC_NP])
{
    int i;

    /* MATLAB interpretation:
     *   E = ones(Np,1)
     */
    for (i = 0; i < MPC_NP; ++i) {
        E[i] = 1.0;
    }
}

static void build_H(double H[MPC_NP][MPC_DELTA_U_SIZE])
{
    int i, j;

    /* MATLAB interpretation for m = 1:
     *   H = tril(ones(Np))
     */
    for (i = 0; i < MPC_NP; ++i) {
        for (j = 0; j < MPC_DELTA_U_SIZE; ++j) {
            H[i][j] = (j <= i) ? 1.0 : 0.0;
        }
    }
}