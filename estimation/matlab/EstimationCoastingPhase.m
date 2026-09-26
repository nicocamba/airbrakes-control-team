function ctx = EstimationCoastingPhase(x, ctx)
%ESTIMATIONCOASTINGPHASE One Kalman estimation step during coasting phase
%
% INPUTS:
%   x   -> current real state [altitude; velocity]
%   ctx -> controller/estimator context
%
% OUTPUT:
%   ctx -> updated context with new estimated state

    %% Format input
    x = x(:);

    %% Current index
    s = ctx.s;

    % Safety check
    if s > ctx.steps
        return;
    end

    %% Read required variables from ctx
    kalman_filter = ctx.kalman_filter;

    var_pressure = ctx.var_pressure;
    gravity = ctx.gravity;
    var_speed = ctx.var_speed;
    bias_var = ctx.bias_var;
    var_bias0 = ctx.var_bias0;
    h = ctx.h;
    var_IMU = ctx.var_IMU;
    sAcc = ctx.sAcc;
    xprev = ctx.xprev;

    %% Compute current air density for barometer variance
    % Here we use the real altitude x(1). If you prefer, this can be based on the
    % previous estimated altitude instead.
    air_dens = interp1( ...
        ctx.air_dens_alt(2,:), ...
        ctx.air_dens_alt(1,:), ...
        x(1), ...
        'linear', ...
        'extrap');

    %% Update barometer measurement variance
    var_baro = var_pressure / (air_dens * gravity)^2;
    kalman_filter.var_baro = var_baro;
    kalman_filter.R = diag([var_baro, var_speed]);

    %% IMU bias random walk
    bias_var = bias_var + randn * sqrt(var_bias0);

    %% Simulated IMU input
    % Measured vertical acceleration from finite difference + noise + bias
    u_Kal = (x(2) - xprev(2)) / h + randn * sqrt(var_IMU) + bias_var;

    %% Prediction step
    kalman_filter.prediction_step(u_Kal);

    %% Simulated measurement vector
    y_Kal = [ ...
        x(1) + randn * sqrt(var_baro); ...
        x(2) + randn * sAcc ...
    ];

    %% Update step
    kalman_filter.update_step(y_Kal);

    %% Save estimated state
    % Depending on your class, x_bar may be the predicted state and x_hat the updated one.
    % Keep x_bar if that is what your implementation uses as the corrected estimate.
    y_est_now = kalman_filter.x_bar;

    if s + 1 <= size(ctx.y_estim, 2)
        ctx.y_estim(:, s+1) = y_est_now;
    end

    %% Save useful estimator data
    ctx.kalman_filter = kalman_filter;
    ctx.bias_var = bias_var;
    ctx.var_baro = var_baro;
    ctx.last_air_dens_est = air_dens;
    ctx.last_u_Kal = u_Kal;
    ctx.last_y_Kal = y_Kal;
    ctx.y_est = y_est_now;
end