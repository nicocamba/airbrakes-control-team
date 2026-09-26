# MATLAB estimator and controller

This directory contains the current MATLAB coasting-phase state estimator/controller and its helper functions. `kalmanCoastingPhaseInit.m` creates controller state, `EstimationCoastingPhase.m` updates the filter, and `kalmanCoastingPhaseStep.m` connects the estimator and controller for a simulation step.

The RocketPy closed-loop notebook starts MATLAB Engine, adds this directory to the MATLAB path, and sets it as the working directory before calling the controller. The MATLAB functions load `flight_data.mat` and `thrust_data.mat` by relative filename; keep these datasets beside the `.m` files.

`launch.csv` is the supplied source table. No conversion script for rebuilding `flight_data.mat` is present in the source archive.
