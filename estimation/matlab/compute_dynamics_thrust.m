function [dXdt] = compute_dynamics_thrust(X, air_dens, CD, T, M, noise)
    % States
    altitude = X(1); 
    speed = X(2);

    % Equation of rocket dynamics
    acc = T / M - 9.77 - 0.5 * air_dens * 9.5e-3 * speed^2 / M * CD + noise;

    dXdt = [speed; acc];
end
