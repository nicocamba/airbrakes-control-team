function [dXdt] = compute_dynamics(X, U, air_dens, CD, CDab, mass, Sab, Sref)
    % States
    altitude = X(1); 
    speed = X(2);

    deflection = U;

    % Equation of rocket dynamics
    acc = -9.77 - 0.5 * air_dens * Sref * speed^2 / mass * (CD + CDab * Sab/Sref * deflection);

    dXdt = [speed; acc];
end

