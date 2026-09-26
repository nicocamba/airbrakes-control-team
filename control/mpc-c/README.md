# C MPC subsystem

This is a modular C port of the coasting-phase MPC structure. `mpc_controller_step(...)` is its top-level entry point. The caller supplies estimated state, references, density/drag tables, previous control, model parameters, weights and input bounds. This code does not read sensors or operate hardware.

Modules:

- `mpc_model`: model parameters, linearization and dynamics.
- `mpc_prediction`: prediction/reference matrices.
- `mpc_constraints`: linear input constraints.
- `mpc_cost`: quadratic-program cost matrices.
- `mpc_solver`: current feasible-warm-start/zero fallback placeholder.
- `mpc_controller`: one-step orchestration.

**Status:** this is not yet a complete optimizer. `solve_qp(...)` does not minimize the objective, and numerical matching against MATLAB remains pending. The original FreeRTOS avionics sketch was excluded because it is an incomplete, unintegrated example and is not part of this C MPC module.
