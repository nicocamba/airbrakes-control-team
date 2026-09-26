# Sunspear Airbrakes: Control and Simulation Algorithms

Research code for the Sunspear supersonic rocket project. This repository brings together RocketPy flight simulations, MATLAB state-estimation/control code, and an early modular C port of the MPC controller.

> **Project status:** research prototype. The C MPC solver is still a placeholder and has not been numerically validated against MATLAB. These files are not a flight-ready control system.

## What is here

- **RocketPy simulations** in [`simulation/rocketpy/`](simulation/rocketpy/): a MATLAB-backed airbrakes simulation and a fixed-deployment reference run. Both use the same rocket, motor, aerodynamic and launch-environment parameters.
- **MATLAB estimator and controller** in [`estimation/matlab/`](estimation/matlab/): the current coasting-phase Kalman/MPC pipeline, helper functions, and its supplied flight datasets. The closed-loop notebook calls this implementation through MATLAB Engine.
- **C MPC port** in [`control/mpc-c/`](control/mpc-c/): a modular implementation of the model, prediction, constraints, quadratic cost, solver interface and top-level control step. The C folder is an algorithm module, not a complete avionics application.

The simulation notebooks call the MATLAB implementation; they do **not** call the C port. The two implementations should not be treated as numerically equivalent.

## Repository layout

```text
control/mpc-c/                 Modular C MPC subsystem
estimation/matlab/             Current MATLAB estimator/controller and input data
simulation/rocketpy/           RocketPy notebooks and shared simulation assets
  assets/                       Motor thrust curve and aerodynamic drag data
  closed_loop_matlab.ipynb      RocketPy simulation using MATLAB control
  open_loop_reference.ipynb     RocketPy reference with fixed airbrake deployment
  results/                      Created by notebook runs; ignored by Git
requirements.txt                Python dependency for the RocketPy notebooks
```

## Running the RocketPy notebooks

Use Python 3.10 or newer. From the repository root, create an environment and install the pinned RocketPy version:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install jupyterlab
jupyter lab
```

Run `simulation/rocketpy/open_loop_reference.ipynb` for the fixed-input reference. It saves CSV/KML outputs under `simulation/rocketpy/results/`.

`simulation/rocketpy/closed_loop_matlab.ipynb` additionally requires a MATLAB installation with MATLAB Engine for Python installed for a compatible Python/MATLAB version. MATLAB Engine is not included in `requirements.txt`; follow the installation instructions for the MATLAB release installed on your machine. The notebook locates the MATLAB source and assets relative to the repository, so it has no machine-specific absolute paths.

The notebooks use a launch-site coordinate and simulation parameters from the source project. Results may depend on RocketPy, MATLAB, and atmospheric-model versions.

## MATLAB data

The controller loads `flight_data.mat` and uses the supplied lookup tables and thrust data. The corresponding CSV source data is retained. The source archive did not include a documented script that regenerates these MAT files from CSV, so the MAT files are included as inputs rather than build products.

## C MPC status

The C modules separate model linearization, prediction matrices, constraints, cost construction, solver and controller orchestration. The solver currently selects a feasible warm start or zero fallback; it does not optimize the quadratic objective. The source project's C notes also identify MATLAB-vs-C numerical validation as outstanding. Integrate and validate the solver before treating its output as an optimized MPC command.

## Provenance and publication

The source archive identifies this material as property of the Sunspear supersonic rocket project. That attribution is preserved here. Confirm that you have permission to publish the team's code and datasets before making this repository public. No open-source license is included because the source material does not grant one.

