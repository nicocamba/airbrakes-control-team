function [dXdt] = compute_dynamics_pert(X, U, air_dens, CD, CDab, noise)
    % States
    altitude = X(1); 
    speed = X(2);

    deflection = U;

    % Equation of rocket dynamics
    acc = -9.77 - 0.5 * air_dens * 9.5e-3 * speed^2 / 16.4813 * (CD + CDab * 4e-3/9.5e-3 * deflection) + noise;

    dXdt = [speed; acc];
end

