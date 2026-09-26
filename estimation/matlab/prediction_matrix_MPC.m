function [Phi_extend, Gamma_extend, Phiref, REF] = prediction_matrix_MPC(n, m, Np, A_c, Phi, x, Gamma_hat, Gamma, Fref)
    % Extended to -->Np
    Phi_extend = zeros(Np * size(A_c,1), size(A_c,2));
    for i = 1:Np
        Phi_extend((i-1)*n+1:i*n, :) = (Phi^i);
    end
    
    Gamma_extend = zeros(Np * n, Np * m);
    
    for i = 1:Np
        for j = 1:i
            sum_phi = zeros(n);  % Suma de potencias de Phi
            for k = 0:(i-j)
                sum_phi = sum_phi + Phi^k;
            end
            row_idx = (i-1)*n + 1 : i*n;
            col_idx = (j-1)*m + 1 : j*m;
            Gamma_extend(row_idx, col_idx) = sum_phi * Gamma;
        end
    end


    Phiref = zeros(n*Np, n);
    Phi_sum = zeros(n);
    for i = 1:Np
        row_idx = (i-1)*n + 1 : i*n;
        col_idx = 1 : n;
        Phi_sum = Phi_sum + Phi^(i-1);
        Phiref(row_idx, col_idx) = Phi_sum;
    end
    
    REF = (eye(n) - Phi) * x + Gamma_hat * Fref;
end