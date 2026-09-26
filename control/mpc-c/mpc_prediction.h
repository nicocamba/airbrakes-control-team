#ifndef MPC_PREDICTION_H
#define MPC_PREDICTION_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mpc_types.h"

/* ============================================================
 * mpc_prediction.h
 *
 * Prediction matrices for MPC.
 *
 * MATLAB tracking:
 *   [Phi_extend, Gamma_extend, Phiref, REF] =
 *       prediction_matrix_MPC(n, m, Np, A_c, Phi, x,
 *                             Gamma_hat, Gamma, Fref);
 * ============================================================ */

int prediction_matrix_MPC(
    const MPC_StateEstimate *y_estim,
    const double A_c[2][2],
    const double Phi[2][2],
    const double Gamma_hat[2][2],
    const double Gamma[2],
    const double Fref[2],
    double Phi_extend[MPC_RP_SIZE][MPC_N_STATES],
    double Gamma_extend[MPC_RP_SIZE][MPC_DELTA_U_SIZE],
    double Phiref[MPC_RP_SIZE][MPC_N_STATES],
    double REF[MPC_N_STATES]
);

#ifdef __cplusplus
}
#endif

#endif /* MPC_PREDICTION_H */