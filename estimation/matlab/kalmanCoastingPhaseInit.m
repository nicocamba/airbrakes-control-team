function ctx = kalmanCoastingPhaseInit(x0)
%KALMANCOASTINGPHASEINIT Initialize context for coasting-phase estimation + MPC
%
% INPUT:
%   x0  -> initial real state at controller start
%          x0 = [altitude; vertical_velocity]
%
% OUTPUT:
%   ctx -> structure containing all persistent data needed by the step function

    %% Input checks
    if nargin < 1
        error('kalmanCoastingPhaseInit requires x0 = [altitude; velocity].');
    end

    if numel(x0) ~= 2
        error('x0 must be a 2x1 or 1x2 vector: [altitude; velocity].');
    end

    x0 = x0(:);  % force column vector

    %% READING AND LOADING FLIGHT DATA
    % From this table of values, we will take in flight parameters
    % (density, drag coefficient) to avoid online computing
    flight_structure = load("flight_data.mat");
    flight_data = flight_structure.flight_data;

    [~, apog_pos] = max(flight_data.Z_m);

    % CONTROL INITIAL INSTANT: M<1 --> Secure zone M<0.95
    start_pos = find(flight_data.Time_s == 16.6, 1, 'first');

    if isempty(start_pos)
        error('No valid start_pos found: check flight_data.Mach and flight_data.Time_s.');
    end

    %% Physical constants
    gravity = 9.77;
    S = 9.5e-3;     % Sref of the rocket

    %% AIRBRAKES DATA
    Sab = 1e-1;
    ab_drag_coeff = 0.25;

    %% Filtered data for density and drag coefficient lookup
    % Stored as 2xN:
    % row 1 -> variable value
    % row 2 -> lookup axis
    air_dens_alt = [ ...
        flight_data.Density_alt(start_pos:apog_pos)'; ...
        flight_data.Z_m(start_pos:apog_pos)' ...
    ];

    CD_speed = [ ...
        flight_data.Total_speed_m_s(start_pos:apog_pos)'; ...
        flight_data.CD_flight(start_pos:apog_pos)' ...
    ];

    %% Linear system characteristics
    n = 2;  % State variables
    m = 1;  % Input variables

    %% Controller parameters
    Np = 10;     % Predictive horizon
    h  = 0.1;    % Sample time

    time = 0:h:(flight_data.Time_s(apog_pos) - flight_data.Time_s(start_pos));
    steps = length(time);

    %% Cost matrices
    Q_step = [1/10^3, 1/10^2];   % States weight
    Q = diag(repmat(Q_step, 1, Np));

    R = 1;

    %% Constraints
    U_min = 0;
    U_max = 1;

    Umin = zeros(Np, 1);
    Umax = ones(Np, 1);

    %% Noise generator for perturbations
    processNoiseStd = 0.2;

    %% Intermediate nominal trajectory --> u = 0.5 constant
    % Assumes original flight_data sampling is about 0.01 s, while controller h = 0.1 s
    steps_num = length(time) + Np;
    altitude_ref = zeros(1, steps_num);
    velocity_ref = zeros(1, steps_num);

    stride = 10;  % 0.1 / 0.01 = 10, keep your original assumption

    for k = 1:steps_num
        idx = start_pos - 1 + stride * k;

        % Saturate index to avoid exceeding table size
        if idx > height(flight_data)
            idx = height(flight_data);
        end

        altitude_ref(k) = flight_data.Z_m(idx);
        velocity_ref(k) = flight_data.Z_speed_m_s(idx);
    end

    [z_apog_nom, z_apog_ind] = max(altitude_ref);

    for k = 1:length(altitude_ref)
        if velocity_ref(k) < 0
            altitude_ref(k) = z_apog_nom;
            velocity_ref(k) = 0;
        end
    end

    %% INITIAL VARIABLES
    % State vector and control value
    X0 = x0;
    U0 = 0;

    x = X0;
    U = [U0, zeros(m, steps-1)];
    u = U0;

    delta_U = zeros(m*Np, 1);

    %% Saving simulation variables
    AB_deployment_sol = zeros(1, steps);
    AB_deployment_sol(1) = U0;

    position_sol = zeros(1, steps+1);
    position_sol(1) = X0(1);

    velocity_sol = zeros(1, steps+1);
    velocity_sol(1) = X0(2);

    %% Optimization timing
    tiempo_total = 0;
    tiempo_iter_i = zeros(steps, 1);

    %% Kalman filter initialization
    % Initial estimated response (adding some error)
    h_0 = x(1) - 20;
    v_0 = x(2) - 2;

    %% Sensors characteristics

    % IMU
    BW = 10;
    dens_noise_IMU = 230; % micro g/sqrt(Hz)
    var_IMU = (dens_noise_IMU * 1e-6 * gravity * sqrt(BW))^2;

    % BIAS
    var_bias  = (50e-03 * gravity / 3)^2;
    var_bias0 = 1e-4;
    bias_var  = 0.1;

    % GPS
    sAcc = 0.1;    % m/s
    var_speed = sAcc^2;

    % Barometer
    accuracyPressure = 250; % Provided by avionics
    var_pressure = (accuracyPressure / 3)^2; % assuming it is a 3*sigma value
    lambda = -6.5e-3;
    Rgas   = 287;
    T0     = 288.15;
    rho0   = 1.225;

    density = rho0 * (1 + lambda * h_0 / T0)^(-1 - gravity / (lambda * Rgas));
    var_baro = var_pressure / (density * gravity)^2;

    y_estim = zeros(3, steps+1);
    y_estim(:,1) = [h_0; v_0; 0];

    kalman_filter = Kalman_develop;
    kalman_filter.var_model = processNoiseStd^2;
    kalman_filter.var_IMU = var_IMU;
    kalman_filter.var_baro = var_baro;
    kalman_filter.var_bias0 = var_bias0;
    kalman_filter.x_0 = [h_0; v_0; 0];
    kalman_filter.P_0 = [20^2, 0, 0;
                         0, 2^2, 0;
                         0, 0, var_bias];
    kalman_filter.Ts = h;
    kalman_filter.H = [1 0 0;
                       0 1 0];
    kalman_filter.R = diag([var_baro, var_speed]);
    kalman_filter.init_filter;

    %% Linearization parameters
    mass = 16.4813;
    rho_der = -0.07e-3;
    Cd_ab = 0.23;
    r_p = zeros(n*Np, 1);
    r_p_log = zeros(n*Np, steps);
    
    %% Internal step counters / memory
    s = 1;
    xprev = X0;
    uprev = U0;
    delta_Uprev = delta_U;

    %% Build output context
    ctx = struct();

    % Raw data
    ctx.flight_data = flight_data;
    ctx.start_pos = start_pos;
    ctx.apog_pos = apog_pos;

    % Physical parameters
    ctx.gravity = gravity;
    ctx.S = S;
    ctx.Sab = Sab;
    ctx.ab_drag_coeff = ab_drag_coeff;
    ctx.mass = mass;
    ctx.rho_der = rho_der;
    ctx.Cd_ab = Cd_ab;

    % Lookup tables
    ctx.air_dens_alt = air_dens_alt;
    ctx.CD_speed = CD_speed;

    % Dimensions and MPC parameters
    ctx.n = n;
    ctx.m = m;
    ctx.Np = Np;
    ctx.h = h;
    ctx.time = time;
    ctx.steps = steps;
    ctx.Q = Q;
    ctx.R = R;
    ctx.U_min = U_min;
    ctx.U_max = U_max;
    ctx.Umin = Umin;
    ctx.Umax = Umax;
    ctx.r_p = r_p;
    ctx.r_p_log = r_p_log;
    % Noise / sensor parameters
    ctx.processNoiseStd = processNoiseStd;
    ctx.BW = BW;
    ctx.dens_noise_IMU = dens_noise_IMU;
    ctx.var_IMU = var_IMU;
    ctx.var_bias = var_bias;
    ctx.var_bias0 = var_bias0;
    ctx.bias_var = bias_var;
    ctx.sAcc = sAcc;
    ctx.var_speed = var_speed;
    ctx.accuracyPressure = accuracyPressure;
    ctx.var_pressure = var_pressure;
    ctx.lambda = lambda;
    ctx.Rgas = Rgas;
    ctx.T0 = T0;
    ctx.rho0 = rho0;
    ctx.var_baro = var_baro;

    % Reference trajectory
    ctx.altitude_ref = altitude_ref;
    ctx.velocity_ref = velocity_ref;
    ctx.z_apog_nom = z_apog_nom;
    ctx.z_apog_ind = z_apog_ind;
    ctx.stride = stride;

    % States and control memory
    ctx.X0 = X0;
    ctx.U0 = U0;
    ctx.x = x;
    ctx.xprev = xprev;
    ctx.u = u;
    ctx.uprev = uprev;
    ctx.U = U;
    ctx.delta_U = delta_U;
    ctx.delta_Uprev = delta_Uprev;
    ctx.s = s;

    % Estimation memory
    ctx.h_0 = h_0;
    ctx.v_0 = v_0;
    ctx.y_estim = y_estim;
    ctx.kalman_filter = kalman_filter;

    % Logging
    ctx.AB_deployment_sol = AB_deployment_sol;
    ctx.position_sol = position_sol;
    ctx.velocity_sol = velocity_sol;
    ctx.tiempo_total = tiempo_total;
    ctx.tiempo_iter_i = tiempo_iter_i;

    % Activation / mode flags
    ctx.controller_initialized = true;
    ctx.controller_enabled = false;  % enable later in step when t>10 && mach<0.95
end