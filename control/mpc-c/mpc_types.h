#ifndef MPC_TYPES_H
#define MPC_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

#define MPC_N_STATES   2
#define MPC_N_INPUTS   1
#define MPC_NP         10

#define MPC_RP_SIZE        (MPC_N_STATES * MPC_NP)
#define MPC_DELTA_U_SIZE   (MPC_N_INPUTS * MPC_NP)

typedef struct
{
    double z;
    double v;
} MPC_StateEstimate;

typedef struct
{
    double gravity;
    double S;
    double Sab;
    double ab_drag_coeff;
    double mass;
    double rho_der;
    double h;
} MPC_Params;

typedef struct
{
    double q_z;
    double q_v;
    double r_du;
} MPC_Weights;

typedef struct
{
    double U_min;
    double U_max;
} MPC_Bounds;

typedef struct
{
    double u;
    double delta_U[MPC_DELTA_U_SIZE];
    int solver_status;
} MPC_Output;

#ifdef __cplusplus
}
#endif

#endif /* MPC_TYPES_H */