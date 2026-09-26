%% Integral correspondiente a la solución particular de la integración de la respuesta del sistema discretizado
function integral_gamma = gamma_integral(A, h)
    n = size(A,1);
    M = [A, eye(n);
         zeros(n), zeros(n)];
    EM = expm(M * h);
    integral_gamma = EM(1:n, n+1:end);  % Esta es ∫₀ʰ e^{Aτ} dτ
end
