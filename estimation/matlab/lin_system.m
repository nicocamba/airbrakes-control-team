function [A_c, B_c] = lin_system(Cd_global, rho_der, S, mass, air_dens, x, Sab, Cd_ab)

A_c = [0, 1;
    -0.5/mass*rho_der*x(2)^2*S*Cd_global, -1/mass*air_dens*x(2)*S*Cd_global];
B_c = [0; -0.5/mass*air_dens*x(2)^2*Sab*Cd_ab];
end