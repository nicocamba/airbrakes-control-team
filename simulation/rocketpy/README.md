# RocketPy simulations

The notebooks use shared inputs from `assets/` and write generated files into `results/`. Launch Jupyter from either this directory or the repository root; the notebooks resolve these paths automatically.

- `closed_loop_matlab.ipynb`: simulates airbrakes controlled by the MATLAB estimator/MPC via MATLAB Engine.
- `open_loop_reference.ipynb`: comparison run with a fixed airbrake deployment level of 0.5.

Install Python dependencies from the repository root using `requirements.txt`. The first notebook also needs a compatible local MATLAB Engine installation.
