#ifndef MPC_SOLVER_H
#define MPC_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mpc_types.h"

/* ============================================================
 * mpc_solver.h
 *
 * QP solver abstraction for MPC.
 *
 * MATLAB tracking:
 *   [delta_U,...] = fmincon(...)
 *
 * Current C implementation:
 * - placeholder / scaffold solver
 * - keeps architecture complete
 * - to be replaced later by a real QP solver
 *   (e.g. OSQP, qpOASES, etc.)
 * ============================================================ */

/* Return codes:
 *   0  -> success
 *  -1  -> invalid input
 *   1  -> warm start accepted as feasible placeholder solution
 *   2  -> warm start infeasible, zero vector used as fallback
 */
int solve_qp(
    const double P[MPC_DELTA_U_SIZE][MPC_DELTA_U_SIZE],
    const double q[MPC_DELTA_U_SIZE],
    const double A_cons[2 * MPC_NP][MPC_DELTA_U_SIZE],
    const double b_cons[2 * MPC_NP],
    const double delta_U0[MPC_DELTA_U_SIZE],
    double delta_U_opt[MPC_DELTA_U_SIZE]
);

#ifdef __cplusplus
}
#endif

#endif /* MPC_SOLVER_H */