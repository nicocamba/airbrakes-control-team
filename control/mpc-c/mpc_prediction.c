#include "mpc_prediction.h"

#include <stddef.h>

/* ============================================================
 * mpc_prediction.c
 *
 * Prediction matrices for MPC.
 *
 * MATLAB source tracking:
 *   kalmanCoastingPhase.m
 * ============================================================ */

/* ============================================================
 * Private helper functions
 * ============================================================ */

/* ------------------------------------------------------------
 * 2x2 identity matrix
 * ------------------------------------------------------------ */
static void mat2_identity(double A[2][2]);

/* ------------------------------------------------------------
 * 2x2 zero matrix
 * ------------------------------------------------------------ */
static void mat2_zero(double A[2][2]);

/* ------------------------------------------------------------
 * 2x2 matrix copy
 * ------------------------------------------------------------ */
static void mat2_copy(const double A[2][2], double B[2][2]);

/* ------------------------------------------------------------
 * 2x2 matrix multiplication
 * ------------------------------------------------------------ */
static void mat2_mul(const double A[2][2], const double B[2][2], double C[2][2]);

/* ------------------------------------------------------------
 * 2x2 matrix addition in place
 * ------------------------------------------------------------ */
static void mat2_add_inplace(double A[2][2], const double B[2][2]);

/* ------------------------------------------------------------
 * 2x2 matrix power
 *
 * MATLAB tracking:
 *   Phi^i
 * ------------------------------------------------------------ */
static void mat2_pow(const double A[2][2], int p, double Ap[2][2]);

/* ------------------------------------------------------------
 * 2x2 matrix times 2x1 vector
 * ------------------------------------------------------------ */
static void mat2_vec2_mul(const double A[2][2], const double x[2], double y[2]);

/* ------------------------------------------------------------
 * Write 2x2 block into big matrix
 * ------------------------------------------------------------ */
static void write_block_2x2(double M[MPC_RP_SIZE][MPC_DELTA_U_SIZE], int row0, int col0, const double B[2][2]);

/* ------------------------------------------------------------
 * Write 2x1 block into big matrix
 * ------------------------------------------------------------ */
static void write_block_2x1(double M[MPC_RP_SIZE][MPC_DELTA_U_SIZE], int row0, int col0, const double b[2]);

/* ------------------------------------------------------------
 * Write 2x2 block into [MPC_RP_SIZE x MPC_N_STATES]
 * ------------------------------------------------------------ */
static void write_block_2x2_state(double M[MPC_RP_SIZE][MPC_N_STATES], int row0, const double B[2][2]);

/* ============================================================
 * Public functions
 * ============================================================ */

int prediction_matrix_MPC(const MPC_StateEstimate *y_estim, const double A_c[2][2], const double Phi[2][2], const double Gamma_hat[2][2], const double Gamma[2], const double Fref[2], double Phi_extend[MPC_RP_SIZE][MPC_N_STATES], double Gamma_extend[MPC_RP_SIZE][MPC_DELTA_U_SIZE], double Phiref[MPC_RP_SIZE][MPC_N_STATES], double REF[MPC_N_STATES])
{
    double x[2];
    double Phi_pow[2][2];
    double Phi_sum[2][2];
    int i, j, k, r, c;

    if (y_estim == NULL || A_c == NULL || Phi == NULL || Gamma_hat == NULL || Gamma == NULL ||
        Fref == NULL || Phi_extend == NULL || Gamma_extend == NULL || Phiref == NULL || REF == NULL) {
        return -1;
    }

    x[0] = y_estim->z;
    x[1] = y_estim->v;

    /* ========================================================
     * MATLAB block:
     *   Phi_extend = zeros(Np * size(A_c,1), size(A_c,2));
     * ======================================================== */
    for (r = 0; r < MPC_RP_SIZE; ++r) {
        for (c = 0; c < MPC_N_STATES; ++c) {
            Phi_extend[r][c] = 0.0;
            Phiref[r][c] = 0.0;
        }
    }

    for (r = 0; r < MPC_RP_SIZE; ++r) {
        for (c = 0; c < MPC_DELTA_U_SIZE; ++c) {
            Gamma_extend[r][c] = 0.0;
        }
    }

    /* ========================================================
     * MATLAB block:
     *   for i = 1:Np
     *       Phi_extend((i-1)*n+1:i*n, :) = (Phi^i);
     *   end
     * ======================================================== */
    for (i = 1; i <= MPC_NP; ++i) {
        int row0 = (i - 1) * MPC_N_STATES;
        mat2_pow(Phi, i, Phi_pow);
        write_block_2x2_state(Phi_extend, row0, Phi_pow);
    }

    /* ========================================================
     * MATLAB block:
     *   Gamma_extend = zeros(Np * n, Np * m);
     *
     *   for i = 1:Np
     *       for j = 1:i
     *           sum_phi = zeros(n);
     *           for k = 0:(i-j)
     *               sum_phi = sum_phi + Phi^k;
     *           end
     *           row_idx = (i-1)*n + 1 : i*n;
     *           col_idx = (j-1)*m + 1 : j*m;
     *           Gamma_extend(row_idx, col_idx) = sum_phi * Gamma;
     *       end
     *   end
     * ======================================================== */
    for (i = 1; i <= MPC_NP; ++i) {
        for (j = 1; j <= i; ++j) {
            double sum_phi[2][2];
            double block[2];

            mat2_zero(sum_phi);

            for (k = 0; k <= (i - j); ++k) {
                mat2_pow(Phi, k, Phi_pow);
                mat2_add_inplace(sum_phi, Phi_pow);
            }

            mat2_vec2_mul(sum_phi, Gamma, block);
            write_block_2x1(Gamma_extend, (i - 1) * MPC_N_STATES, (j - 1) * MPC_N_INPUTS, block);
        }
    }

    /* ========================================================
     * MATLAB block:
     *   Phiref = zeros(n*Np, n);
     *   Phi_sum = zeros(n);
     *   for i = 1:Np
     *       row_idx = (i-1)*n + 1 : i*n;
     *       Phi_sum = Phi_sum + Phi^(i-1);
     *       Phiref(row_idx, :) = Phi_sum;
     *   end
     * ======================================================== */
    mat2_zero(Phi_sum);

    for (i = 1; i <= MPC_NP; ++i) {
        int row0 = (i - 1) * MPC_N_STATES;

        mat2_pow(Phi, i - 1, Phi_pow);
        mat2_add_inplace(Phi_sum, Phi_pow);
        write_block_2x2_state(Phiref, row0, Phi_sum);
    }

    /* ========================================================
     * MATLAB block:
     *   REF = (eye(n) - Phi) * x + Gamma_hat * Fref;
     * ======================================================== */
    {
        double I[2][2];
        double I_minus_Phi[2][2];
        double term1[2];
        double term2[2];

        mat2_identity(I);

        for (r = 0; r < 2; ++r) {
            for (c = 0; c < 2; ++c) {
                I_minus_Phi[r][c] = I[r][c] - Phi[r][c];
            }
        }

        mat2_vec2_mul(I_minus_Phi, x, term1);
        mat2_vec2_mul(Gamma_hat, Fref, term2);

        REF[0] = term1[0] + term2[0];
        REF[1] = term1[1] + term2[1];

        /* MATLAB REF is 2x1, but in your controller scaffold REF was
         * sized as MPC_RP_SIZE. To keep compatibility with that scaffold,
         * we fill only the first 2 positions and zero the rest.
         */
      
        
    }

    return 0;
}

/* ============================================================
 * Private helper implementations
 * ============================================================ */

static void mat2_identity(double A[2][2])
{
    A[0][0] = 1.0; A[0][1] = 0.0;
    A[1][0] = 0.0; A[1][1] = 1.0;
}

static void mat2_zero(double A[2][2])
{
    A[0][0] = 0.0; A[0][1] = 0.0;
    A[1][0] = 0.0; A[1][1] = 0.0;
}

static void mat2_copy(const double A[2][2], double B[2][2])
{
    B[0][0] = A[0][0]; B[0][1] = A[0][1];
    B[1][0] = A[1][0]; B[1][1] = A[1][1];
}

static void mat2_mul(const double A[2][2], const double B[2][2], double C[2][2])
{
    double T[2][2];

    T[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
    T[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
    T[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
    T[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];

    mat2_copy(T, C);
}

static void mat2_add_inplace(double A[2][2], const double B[2][2])
{
    A[0][0] += B[0][0]; A[0][1] += B[0][1];
    A[1][0] += B[1][0]; A[1][1] += B[1][1];
}

static void mat2_pow(const double A[2][2], int p, double Ap[2][2])
{
    int i;
    double temp[2][2];
    double result[2][2];

    if (p == 0) {
        mat2_identity(Ap);
        return;
    }

    mat2_identity(result);

    for (i = 0; i < p; ++i) {
        mat2_mul(result, A, temp);
        mat2_copy(temp, result);
    }

    mat2_copy(result, Ap);
}

static void mat2_vec2_mul(const double A[2][2], const double x[2], double y[2])
{
    y[0] = A[0][0] * x[0] + A[0][1] * x[1];
    y[1] = A[1][0] * x[0] + A[1][1] * x[1];
}

static void write_block_2x2(double M[MPC_RP_SIZE][MPC_DELTA_U_SIZE], int row0, int col0, const double B[2][2])
{
    M[row0 + 0][col0 + 0] = B[0][0];
    M[row0 + 0][col0 + 1] = B[0][1];
    M[row0 + 1][col0 + 0] = B[1][0];
    M[row0 + 1][col0 + 1] = B[1][1];
}

static void write_block_2x1(double M[MPC_RP_SIZE][MPC_DELTA_U_SIZE], int row0, int col0, const double b[2])
{
    M[row0 + 0][col0] = b[0];
    M[row0 + 1][col0] = b[1];
}

static void write_block_2x2_state(double M[MPC_RP_SIZE][MPC_N_STATES], int row0, const double B[2][2])
{
    M[row0 + 0][0] = B[0][0];
    M[row0 + 0][1] = B[0][1];
    M[row0 + 1][0] = B[1][0];
    M[row0 + 1][1] = B[1][1];
}