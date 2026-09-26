#include "mpc_solver.h"

#include <stddef.h>

/* ============================================================
 * mpc_solver.c
 *
 * Placeholder QP solver for MPC.
 *
 * MATLAB source tracking:
 *   kalmanCoastingPhase.m
 *
 * MATLAB block:
 *   [delta_U,...] = fmincon(fun, delta_U0, A, b, [], [], [], [], [], options);
 *
 * Current implementation:
 * - this is NOT a real optimizer yet
 * - it provides a stable architectural placeholder
 * - later it should be replaced by a real QP solver
 * ============================================================ */

/* ============================================================
 * Private helper functions
 * ============================================================ */

/* ------------------------------------------------------------
 * Copy vector
 * ------------------------------------------------------------ */
static void vec_copy(const double *src, double *dst, int n);

/* ------------------------------------------------------------
 * Set vector to zero
 * ------------------------------------------------------------ */
static void vec_zero(double *x, int n);

/* ------------------------------------------------------------
 * Check linear inequalities:
 *
 *   A_cons * x <= b_cons
 * ------------------------------------------------------------ */
static int check_inequalities(const double A_cons[2 * MPC_NP][MPC_DELTA_U_SIZE], const double b_cons[2 * MPC_NP], const double x[MPC_DELTA_U_SIZE]);

/* ------------------------------------------------------------
 * Evaluate quadratic cost:
 *
 *   J(x) = 1/2 x' P x + q' x
 *
 * Not strictly needed for the current placeholder, but useful
 * for debugging and future extension.
 * ------------------------------------------------------------ */
static double evaluate_qp_cost(const double P[MPC_DELTA_U_SIZE][MPC_DELTA_U_SIZE], const double q[MPC_DELTA_U_SIZE], const double x[MPC_DELTA_U_SIZE]);

/* ============================================================
 * Public functions
 * ============================================================ */

int solve_qp(const double P[MPC_DELTA_U_SIZE][MPC_DELTA_U_SIZE], const double q[MPC_DELTA_U_SIZE], const double A_cons[2 * MPC_NP][MPC_DELTA_U_SIZE], const double b_cons[2 * MPC_NP], const double delta_U0[MPC_DELTA_U_SIZE], double delta_U_opt[MPC_DELTA_U_SIZE])
{
    double zero_candidate[MPC_DELTA_U_SIZE];
    int feasible_warm;
    int feasible_zero;

    if (P == NULL || q == NULL || A_cons == NULL || b_cons == NULL || delta_U_opt == NULL) {
        return -1;
    }

    /* ========================================================
     * MATLAB tracking:
     *   fmincon starts from delta_U0 as initial guess
     *
     * Current placeholder logic:
     * 1. If warm start is feasible, use it directly
     * 2. Else, try zero vector
     * 3. Else, still return zero vector as emergency fallback
     * ======================================================== */

    vec_zero(zero_candidate, MPC_DELTA_U_SIZE);

    feasible_warm = 0;
    feasible_zero = 0;

    if (delta_U0 != NULL) {
        feasible_warm = check_inequalities(A_cons, b_cons, delta_U0);
    }

    feasible_zero = check_inequalities(A_cons, b_cons, zero_candidate);

    if (delta_U0 != NULL && feasible_warm) {
        vec_copy(delta_U0, delta_U_opt, MPC_DELTA_U_SIZE);

        /* Optional debug hook:
         * evaluate_qp_cost(P, q, delta_U_opt);
         */
        return 1;
    }

    if (feasible_zero) {
        vec_copy(zero_candidate, delta_U_opt, MPC_DELTA_U_SIZE);

        /* Optional debug hook:
         * evaluate_qp_cost(P, q, delta_U_opt);
         */
        return 2;
    }

    /* Emergency fallback:
     * return zero vector even if constraints look inconsistent.
     * Higher layer may decide to keep previous control.
     */
    vec_zero(delta_U_opt, MPC_DELTA_U_SIZE);
    return 2;
}

/* ============================================================
 * Private helper implementations
 * ============================================================ */

static void vec_copy(const double *src, double *dst, int n)
{
    int i;

    for (i = 0; i < n; ++i) {
        dst[i] = src[i];
    }
}

static void vec_zero(double *x, int n)
{
    int i;

    for (i = 0; i < n; ++i) {
        x[i] = 0.0;
    }
}

static int check_inequalities(const double A_cons[2 * MPC_NP][MPC_DELTA_U_SIZE], const double b_cons[2 * MPC_NP], const double x[MPC_DELTA_U_SIZE])
{
    int i, j;
    const double tol = 1e-9;

    for (i = 0; i < 2 * MPC_NP; ++i) {
        double lhs = 0.0;

        for (j = 0; j < MPC_DELTA_U_SIZE; ++j) {
            lhs += A_cons[i][j] * x[j];
        }

        if (lhs > b_cons[i] + tol) {
            return 0;
        }
    }

    return 1;
}

static double evaluate_qp_cost(const double P[MPC_DELTA_U_SIZE][MPC_DELTA_U_SIZE], const double q[MPC_DELTA_U_SIZE], const double x[MPC_DELTA_U_SIZE])
{
    int i, j;
    double quad = 0.0;
    double lin = 0.0;

    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        for (j = 0; j < MPC_DELTA_U_SIZE; ++j) {
            quad += x[i] * P[i][j] * x[j];
        }
    }

    for (i = 0; i < MPC_DELTA_U_SIZE; ++i) {
        lin += q[i] * x[i];
    }

    return 0.5 * quad + lin;
}