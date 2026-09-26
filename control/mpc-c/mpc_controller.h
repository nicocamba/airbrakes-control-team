#ifndef MPC_CONTROLLER_H
#define MPC_CONTROLLER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mpc_types.h"

/* ============================================================
 * mpc_controller.h
 *
 * High-level MPC controller interface.
 *
 * MATLAB tracking:
 *   kalmanCoastingPhase.m
 *
 * This module is the top-level entry point of the CONTROL block.
 * It is not the global program main(), but the main interface
 * of the MPC subsystem.
 * ============================================================ */

/* ------------------------------------------------------------
 * Initialize delta_U vector
 *
 * MATLAB tracking:
 *   delta_U = zeros(m*Np,1);
 * ------------------------------------------------------------ */
void mpc_init_delta_u(double delta_U[MPC_DELTA_U_SIZE]);

/* ------------------------------------------------------------
 * Execute one MPC control step
 *
 * MATLAB high-level tracking:
 *
 *   for s = 1:total_steps
 *       xprev = x;
 *       uprev = u;
 *       delta_Uprev = delta_U;
 *
 *       for i = 1:Np
 *           r_p(n*i-1:n*i) = [altitude_ref(s+i-1); velocity_ref(s+i-1)];
 *       end
 *
 *       [air_dens, drag_coeff] = compute_current_parameters(...)
 *       Cd_global = drag_coeff + Sab / S * ab_drag_coeff * u
 *       [A_c, B_c] = lin_system(...)
 *       Phi = expm(A_c*h)
 *       Gamma_hat = gamma_integral(A_c, h)
 *       Gamma = Gamma_hat * B_c
 *       Fref = compute_dynamics(...)
 *       [Phi_extend, Gamma_extend, Phiref, REF] = prediction_matrix_MPC(...)
 *       [A, b] = MPC_Lconsfun(...)
 *       fun = @(delta_U) MPC_objfun(...)
 *       [delta_U,...] = fmincon(...)
 *       u = uprev + delta_U(1)
 *   end
 *
 * C architecture:
 * - this function performs one control iteration
 * - it is intended to be called from an external loop
 * ------------------------------------------------------------ */
int mpc_controller_step(
    const MPC_StateEstimate *y_estim,
    const double *altitude_ref,
    const double *velocity_ref,
    int ref_length,
    int s,
    const double (*air_dens_alt)[2],
    int air_dens_alt_len,
    const double (*CD_speed)[2],
    int cd_speed_len,
    double uprev,
    const double delta_Uprev[MPC_DELTA_U_SIZE],
    const MPC_Params *params,
    const MPC_Weights *weights,
    const MPC_Bounds *bounds,
    MPC_Output *out
);

#ifdef __cplusplus
}
#endif

#endif /* MPC_CONTROLLER_H */