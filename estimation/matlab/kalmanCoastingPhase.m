function [input, ctx] = kalmanCoastingPhase(t, x, ctx)
%KALMANCOASTINGPHASE One MPC control step during coasting phase
%
% INPUTS:
%   t    -> current time
%   x    -> current state [altitude; velocity]
%   mach -> current Mach number
%   ctx  -> controller context created by kalmanCoastingPhaseInit
%
% OUTPUTS:
%   input -> commanded airbrakes deployment in [0, 1]
%   ctx   -> updated controller context

    %% Input formatting
    x = x(:);

    %% Read current step index
    s = ctx.s;

    % Safety: if horizon/log arrays are exhausted, keep last command
    if s > ctx.steps
        input = ctx.u;
        return;
    end

    ctx.controller_enabled = true;

    %% Short aliases from ctx
    n = ctx.n;
    m = ctx.m;
    Np = ctx.Np;
    h = ctx.h;

    S = ctx.S;
    Sab = ctx.Sab;
    ab_drag_coeff = ctx.ab_drag_coeff;
    mass = ctx.mass;
    rho_der = ctx.rho_der;

    Q = ctx.Q;
    R = ctx.R;
    Umin = ctx.Umin;
    Umax = ctx.Umax;

    air_dens_alt = ctx.air_dens_alt;
    CD_speed = ctx.CD_speed;

    altitude_ref = ctx.altitude_ref;
    velocity_ref = ctx.velocity_ref;

    uprev = ctx.u;
    delta_Uprev = ctx.delta_U;
    y_estim = ctx.y_estim;

    %% Build reference over prediction horizon
    %r_p = zeros(n*Np, 1);
    %for k = 1:Np
    %    idx_ref = s + k - 1;
%
    %    % Saturate index to end of reference vector
    %    if idx_ref > length(altitude_ref)
    %        idx_ref = length(altitude_ref);
    %    end
%
    %    r_p(n*k-1:n*k) = [altitude_ref(idx_ref); velocity_ref(idx_ref)];
    %end

        %% Reference trajectory over the horizon (MPC THEORY) --> This is the "desired" trajectory over the corresponding seconds.
    r_p = zeros(n*Np, 1);
    for i = 1:Np
        r_p(n*i-1:n*i) = [altitude_ref(s+i-1); velocity_ref(s+i-1)];
    end

    %% Determine current density and drag using estimated state
%    [air_dens, drag_coeff] = compute_current_parameters( ...
%        air_dens_alt, CD_speed, y_estim(1,s+1), y_estim(2,s+1));
    [air_dens, drag_coeff] = compute_current_parameters( ...
        air_dens_alt, CD_speed, x(1), x(2));


    %% Linearization around current estimated point
    Cd_global = drag_coeff + Sab / S * ab_drag_coeff * uprev;

    %[A_c, B_c] = lin_system( ...
    %    Cd_global, rho_der, S, mass, air_dens, y_estim(:,s+1), Sab, ab_drag_coeff);

    [A_c, B_c] = lin_system( ...
        Cd_global, rho_der, S, mass, air_dens, x, Sab, ab_drag_coeff);



    %% Discrete system
    Phi = expm(A_c * h);
    Gamma_hat = gamma_integral(A_c, h);
    Gamma = Gamma_hat * B_c;

   % Fref = compute_dynamics(y_estim(:,s+1), uprev, air_dens, drag_coeff, ab_drag_coeff);
    Fref = compute_dynamics(x, uprev, air_dens, drag_coeff, ab_drag_coeff, mass, Sab, S);
    %% Prediction matrices
    [Phi_extend, Gamma_extend, Phiref, REF] = prediction_matrix_MPC( ...
        n, m, Np, A_c, Phi, x, Gamma_hat, Gamma, Fref);

   %% Optimization with constraints (MPC THEORY)
    E = repmat(eye(m), Np, 1);
    H = zeros(Np*m, Np*m);
    for i = 1:Np
        for j = 1:i
            row_idx = (i-1)*m + 1 : i*m;
            col_idx = (j-1)*m + 1 : j*m;
            H(row_idx, col_idx) = eye(m);
        end
    end

    %% Warm start
    delta_U0 = zeros(m*Np,1);
    for i=1:Np-1
         delta_U0(i) = delta_Uprev(i+1);
    end

    %% Optimization problem
    [A, b] = MPC_Lconsfun(Umin, Umax, uprev, E, H);

    %fun = @(delta_U) MPC_objfun( ...
    %    y_estim(1:2,s+1), Phi_extend, Gamma_extend, r_p, Q, R, delta_U, Phiref, REF);

    fun = @(delta_U) MPC_objfun( ...
         x, Phi_extend, Gamma_extend, r_p, Q, R, delta_U, Phiref, REF);


    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'OutputFcn', @outfun_rocket, ...
        'MaxIterations', 100, ...
        'MaxFunctionEvaluations', 1000, ...
        'OptimalityTolerance', 1e-3, ...
        'StepTolerance', 1e-6, ...
        'FunctionTolerance', 1e-6);


    %% Solve MPC
    tic;
    [delta_U, ~, exitflag] = fmincon(fun, delta_U0, A, b, [], [], [], [], [], options);
    tiempo_iter = toc;

 %  %% Fallback if optimizer fails badly
 %  if isempty(delta_U) || exitflag < 0
 %      delta_U = delta_Uprev;
 %  end

    %% Apply first control increment
    u = uprev + delta_U(1);

    input = u;

    %% Save control memory
    if s + 1 <= size(ctx.U, 2)
        ctx.U(:, s+1) = u;
    end

    ctx.u = u;
    ctx.uprev = u;
    ctx.delta_U = delta_U;
    ctx.delta_Uprev = delta_U;

    %% Save logs
    ctx.AB_deployment_sol(s) = u;

    if s + 1 <= length(ctx.position_sol)
        ctx.position_sol(s+1) = x(1);
        ctx.velocity_sol(s+1) = x(2);
    end

    ctx.tiempo_iter_i(s) = tiempo_iter;
    ctx.tiempo_total = ctx.tiempo_total + tiempo_iter;

    %% Save extra current values if useful
    ctx.r_p = r_p;
    ctx.r_p_log(:, s) = r_p;
    ctx.xprev = x;
    ctx.last_air_dens = air_dens;
    ctx.last_drag_coeff = drag_coeff;
    ctx.last_exitflag = exitflag;

    %% Advance step counter
    ctx.s = s + 1;
end