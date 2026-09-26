%% Restricciones lineales del problema de opitmización
function [A, b] = MPC_Lconsfun(Umin, Umax, u, E, H)

    % Inequality constraints in a column vector
    A = [-H;   % Min U
        H];      % Max U
       % -Gamma_extend_select;     % Min Y
       % Gamma_extend_select];     % Max Y


    b = [(-Umin+E*u);
        (Umax-E*u)];
        %(-Ymin + Phi_extend_select*x + Phiref_select*REF);
        %(Ymax  - Phi_extend_select*x - Phiref_select*REF)];

end