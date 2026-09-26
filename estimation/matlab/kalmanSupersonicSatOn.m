%% ALTITUDE LINEAR MPC CONTROLLER. SUNSPEAR ROCKET. FILTERED SIGNAL
clc;
clear;
close all;

%% READING AND LOADING FLIGHT DATA
flight_structure = load("flight_data.mat");
flight_data = flight_structure.flight_data;

% At t=7.22s, thrust off
start_pos = find(flight_data.Time_s > 7.8);
start_pos = start_pos(1);

thrust_off = find(flight_data.Time_s >= 14.5);
thrust_off = thrust_off(1); 

gravity = 9.77; 

S = 9.5e-3;     % Sref of the rocket

%% READING AND LOADING AB's DATA
Sab = 4e-3;
ab_drag_coeff = 0.22;

%% Filtrated data for density and temperature --> Only in control interval
air_dens_alt = [flight_data.Density_alt(start_pos : thrust_off)' 
    flight_data.Z_m(start_pos : thrust_off)']';

CD_speed = [flight_data.Total_speed_m_s(start_pos : thrust_off)' 
    flight_data.CD_flight(start_pos : thrust_off)']';

%% Linear system characteristics
n = 2; % State variables
m = 1; % Input variables
p = 2; % Output measurable variables

%% Definition of controller parameters       
Np=10;   % Predictive horizon
h=0.1;    % Sample time       

time = 0:h:(flight_data.Time_s(thrust_off) - flight_data.Time_s(start_pos));
steps = length(time);

% COST MATRIX 
Q_step = [1/10^3  1/(10^3)];     % States weight
Q = diag(repmat(Q_step, 1, Np));    

R = 5;

% Constraints
% Max and min AB deployment
U_min = 0;
U_max = 1;

Umin = zeros(Np, 1); % Expand the condition to hole horizont
Umax = ones(Np, 1);


%% Trajectory used as nominal --> u = 0
x0ref = [3835.93; 493];
input = 0;
steps_num = length(time)+100;

[altitude_ref, velocity_ref] = compute_ConstInput_apog(x0ref, input, steps_num, air_dens_alt,...
    CD_speed, ab_drag_coeff, h);

[z_apog_nom, z_apog_ind] = max(altitude_ref);

%% INITIAL VARIABLES
% State vector and control value
X0 = x0ref + [-10;2];    % +10 m deviation and +2m/s in speed
U0 = 0;

% Initialization
x = X0;
U =[U0 zeros(m, steps-1)];
u = U0;

delta_U = zeros(m*Np,1);

AB_deployment_sol = [U0(1,1) zeros(1,steps-1)];
position_sol = [X0(1) zeros(1,steps)];
velocity_sol = [X0(2) zeros(1,steps)];

tiempo_total = 0;

tiempo_iter_i = zeros(steps, 1);


%% Initialization of the kalman filter
% We initialize the estimated response (adding some error)
h_0 = x(1) + 20;
v_0 = x(2) + 5;

% Sensors characteristics
%IMU
BW = 10;
dens_noise_IMU = 230; % Technical specifications: micro g/sqrt(Hz)
var_IMU = (dens_noise_IMU*1e-6*gravity*sqrt(BW))^2;

%BIAS
var_bias = (50e-03*gravity/3)^2;
var_bias0 = 1e-4;
bias_var=0.1;

% SAT
sAcc = 0.2;     % 0.05 provided in Datasheed for low speeds, will depend on velocity
var_speed = sAcc^2;

y_estim = [[h_0; v_0; 0] zeros(3,steps)]; % Initial deviation for filter

kalman_filter = Kalman_develop;
kalman_filter.var_IMU = var_IMU;
kalman_filter.var_speed = var_speed;
kalman_filter.var_bias0 = var_bias0;
kalman_filter.x_0 = [h_0; v_0; 0];
kalman_filter.P_0 = [20^2, 0, 0; 
                    0, 5^2, 0;
                    0, 0, var_bias];      % P0 from Thrust off!
kalman_filter.H = [0 1 0];
kalman_filter.R = var_speed;
kalman_filter.Ts = h;

kalman_filter.init_filter;

%% Linealization parameters
mass = 16.4813;
Cd_ab = 0.23;
r_p = zeros(n*Np,1);

% Noise generator for perturbations
processNoiseStd = 0.2;
total_steps = steps;
noise_seq = processNoiseStd * randn(total_steps, 1);
air_dens = 78181; 

%% Development
for s=1:total_steps
    xprev = x;
    %% Reference trajectory over the horizon (MPC THEORY)
    for i = 1:Np
        r_p(n*i-1:n*i) = [altitude_ref(s+i-1); velocity_ref(s+i-1)];
    end
    %% Determination CD and rho --> Using y of estimation!
    airDensOld = air_dens; 

    [air_dens, drag_coeff] = compute_current_parameters(air_dens_alt, CD_speed,...
        y_estim(1,s), y_estim(2,s));   
    %% Linealization around initial point estimated by kalman

    Cd_global = drag_coeff;
    if s == 1
        rho_der = -9e-5; 
    else 
        heightDiff = y_estim(1,s) - y_estim(1,s-1); 
        rho_der = (air_dens - airDensOld)/heightDiff;
    end 

    [A_c, B_c] = lin_system(Cd_global, rho_der, S, mass, air_dens, y_estim(:,s), Sab, ab_drag_coeff);

    %% CONVERSION TO DISCRETE SYSTEM (MPC THEORY)
    Phi = expm(A_c*h);
    Gamma_hat = gamma_integral(A_c, h);
    
    Gamma = Gamma_hat*B_c;

    Fref = compute_dynamics(y_estim(:,s), u, air_dens, drag_coeff, ab_drag_coeff);

    %% Assembly of prediction matrix over the horizon (MPC THEORY)
    % Extended to the dynamics and the state-->Np
    [Phi_extend, Gamma_extend, Phiref, REF] = prediction_matrix_MPC(n, m, Np, A_c, Phi, x, Gamma_hat, Gamma, Fref);

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

   % Initialization of the optimizer
   delta_U0 = zeros(m*Np,1);       

   for i=1:Np-1
        delta_U0(i) = delta_U(i+1);
   end

   % Optimization problem computation
   [A, b] = MPC_Lconsfun(Umin, Umax, u, E, H);   
   fun = @(delta_U) MPC_objfun(y_estim(1:2,s), Phi_extend, Gamma_extend, r_p, Q, R, delta_U, Phiref, REF);

   options = optimoptions('fmincon','Display', 'iter', 'Algorithm', 'sqp','OutputFcn', @outfun_rocket,...
       'MaxIterations', 100,'MaxFunctionEvaluations', 1000,'OptimalityTolerance', 1e-3,'StepTolerance', 1e-6,'FunctionTolerance', 1e-6);

   tic;
   [delta_U,fval,exitflag,output,lambda,grad,hessian] =  fmincon(fun,delta_U0, A, b, [], [], [], [], [],options);
   tiempo_iter=toc;
   tiempo_total = tiempo_total + tiempo_iter;
   tiempo_iter_i(s) = tiempo_iter;

   u = u + delta_U(1);

   %% Input computed --> Integration of the system
    % Determination CD and rho for real system

   [air_dens_r, drag_coeff_r] = compute_current_parameters(air_dens_alt, CD_speed,...
      x(1), x(2));

   t_span = h*[s s+1];

   noise = noise_seq(s);
   % SEE THAT U=0, WE DON'T HAVE THE INPUT COMPUTATIONS DELETED JUST IN
   % CASE THEY ARE NEEDED IN THE FUTURE (WARM START)
   [t, Xode] = ode45(@(t, X) compute_dynamics_pert(X, 0, air_dens_r, drag_coeff_r, ab_drag_coeff, noise), t_span, x);   
   
   x = Xode(end,:)';

   U(:,s+1) = u(:,1);

   %% Save real simulation
   AB_deployment_sol(s) = U(s);
   velocity_sol(s+1) = x(2);
   position_sol(s+1) = x(1);

    %% State estimation

    var_speed = sAcc^2;      % changing variance or not
    kalman_filter.var_speed = var_speed;
    kalman_filter.R = var_speed;

    bias_var = bias_var + randn*sqrt(var_bias0);   % random walk (o bias constante)
    u_Kal = (x(2) - xprev(2))/h + randn * sqrt(var_IMU) + bias_var; 

    %u_Kal = thrust_curve(s)/mass_T(s) - gravity - 0.5 * air_dens_r * x(2)^2 * S * drag_coeff_r  / mass_T(s) + ...
    %    noise + randn * sqrt(var_IMU);

    kalman_filter.prediction_step(u_Kal);

    y_Kal = x(2) + randn * sAcc;

    kalman_filter.prediction_step(u_Kal);

    kalman_filter.update_step(y_Kal);
    

    %% Save state estimation
    y_estim(:, s + 1) = kalman_filter.x_bar;

end

error = zeros(1, length(position_sol)); 

for i = 1:length(position_sol)
    error(i) = 100*abs(1 - position_sol(i)/altitude_ref(i)); 
end  
meanError = mean(error); 


%% Extra        -- >Selected blind stage: 2.7s --> 7.8s
%size_flight = fix(size(flight_data.Time_s,1)/10)
%acc_step = zeros(size_flight-1,1);
%for i=2:size_flight
%    acc_step(i) = abs((flight_data.Total_speed_m_s(10*i)-flight_data.Total_speed_m_s(10*(i-1)))/(flight_data.Time_s(10*i)-flight_data.Time_s(10*(i-1))))/gravity;
%end
%[acc_max, acc_max_i] = max(acc_step)
%% Plots

figure;
plot(velocity_sol, position_sol, 'r-', 'LineWidth', 1.5); hold on;
plot(velocity_ref(1:steps+1), altitude_ref(1:steps+1), 'b--', 'LineWidth', 1.5)
xlabel('v');
ylabel('z');
title("z -- mean error of " + meanError + " % ");
legend('z(m)', 'z_n');
grid on;

t = (0:steps) * h;

figure('Color','w','Position',[200 200 900 600]);

% --- Altitude error ---
subplot(2,1,1)
plot(t, y_estim(1,:) - position_sol, ...
     'LineWidth', 1.8, 'Color',[0.85 0.1 0.1]);
grid on
ylabel('Altitude error [m]', 'FontSize', 11)
title('Altitude Estimation Error', 'FontSize', 12)
xlim([t(1) t(end)])

% --- Velocity error ---
subplot(2,1,2)
plot(t, y_estim(2,:) - velocity_sol, ...
     'LineWidth', 1.8, 'Color',[0.1 0.2 0.85]);
grid on
xlabel('Time [s]', 'FontSize', 11)
ylabel('Velocity error [m/s]', 'FontSize', 11)
title('Velocity Estimation Error', 'FontSize', 12)
xlim([t(1) t(end)])

% --- Global title ---
sgtitle('Kalman Filter Estimation Error during Supersonic GNSS allowed Phase, T_s=0.1s', ...
        'FontSize', 14, 'FontWeight','bold')