#ifndef MPC_MODEL_H
#define MPC_MODEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mpc_types.h"

/* ============================================================
 * mpc_model.h
 *
 * Model-related functions used by the controller.
 *
 * MATLAB tracking:
 * These functions correspond to model/aerodynamic blocks used in:
 *   kalmanCoastingPhase.m
 * ============================================================ */

/* ------------------------------------------------------------
 * MATLAB tracking:
 *   [air_dens, drag_coeff] = compute_current_parameters(
 *       air_dens_alt, CD_speed, y_estim(1,s), y_estim(2,s));
 * ------------------------------------------------------------ */
int compute_current_parameters(
    const double (*air_dens_alt)[2],
    int air_dens_alt_len,
    const double (*CD_speed)[2],
    int cd_speed_len,
    double altitude,
    double velocity,
    double *air_dens,
    double *drag_coeff
);

/* ------------------------------------------------------------
 * MATLAB tracking:
 *   Cd_global = drag_coeff + Sab / S * ab_drag_coeff * u;
 *   [A_c, B_c] = lin_system(Cd_global, rho_der, S, mass,
 *                           air_dens, x, Sab, ab_drag_coeff);
 * ------------------------------------------------------------ */
int lin_system(
    double Cd_global,
    double rho_der,
    double S,
    double mass,
    double air_dens,
    const MPC_StateEstimate *y_estim,
    double Sab,
    double ab_drag_coeff,
    double A_c[2][2],
    double B_c[2]
);

/* ------------------------------------------------------------
 * MATLAB tracking:
 *   Gamma_hat = gamma_integral(A_c, h);
 * ------------------------------------------------------------ */
int gamma_integral(
    const double A[2][2],
    double h,
    double Gamma_hat[2][2]
);

/* ------------------------------------------------------------
 * MATLAB tracking:
 *   Fref = compute_dynamics(y_estim(:,s), u, air_dens,
 *                           drag_coeff, ab_drag_coeff);
 * ------------------------------------------------------------ */
int compute_dynamics(
    const MPC_StateEstimate *y_estim,
    double U,
    double air_dens,
    double drag_coeff,
    double ab_drag_coeff,
    double Fref[2]
);

#ifdef __cplusplus
}
#endif

#endif /* MPC_MODEL_H */