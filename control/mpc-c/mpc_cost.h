#ifndef MPC_COST_H
#define MPC_COST_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mpc_types.h"

/* ============================================================
 * mpc_cost.h
 *
 * MPC cost function reformulated as QP:
 *
 * MATLAB tracking:
 *   Y = Phi_extend*x + Gamma_extend*delta_U + Phiref*REF;
 *   MPC_objfun = 1/2*(r_p-Y)'*Q*(r_p-Y) + 1/2*delta_U'*R*delta_U;
 *
 * C/QP reformulation:
 *   build P and q such that
 *
 *       1/2 * delta_U' * P * delta_U + q' * delta_U
 *
 * matches the MATLAB objective up to constant terms.
 * ============================================================ */

int MPC_objfun_to_qp(
    const MPC_StateEstimate *y_estim,
    const double Phi_extend[MPC_RP_SIZE][MPC_N_STATES],
    const double Gamma_extend[MPC_RP_SIZE][MPC_DELTA_U_SIZE],
    const double r_p[MPC_RP_SIZE],
    const double Phiref[MPC_RP_SIZE][MPC_N_STATES],
    const double REF[MPC_N_STATES],
    const MPC_Weights *weights,
    double P[MPC_DELTA_U_SIZE][MPC_DELTA_U_SIZE],
    double q[MPC_DELTA_U_SIZE]
);

#ifdef __cplusplus
}
#endif

#endif /* MPC_COST_H */