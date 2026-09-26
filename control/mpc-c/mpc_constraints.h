#ifndef MPC_CONSTRAINTS_H
#define MPC_CONSTRAINTS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mpc_types.h"

/* ============================================================
 * mpc_constraints.h
 *
 * Linear inequality constraints for MPC.
 *
 * MATLAB tracking:
 *   [A, b] = MPC_Lconsfun(Umin, Umax, u, E, H);
 *
 * In this implementation:
 * - E and H are built internally
 * - constraints are returned in MATLAB-like form:
 *
 *       A_cons * delta_U <= b_cons
 * ============================================================ */

int MPC_Lconsfun(double U_min, double U_max, double uprev, double A_cons[2 * MPC_NP][MPC_DELTA_U_SIZE], double b_cons[2 * MPC_NP]);

#ifdef __cplusplus
}
#endif

#endif /* MPC_CONSTRAINTS_H */