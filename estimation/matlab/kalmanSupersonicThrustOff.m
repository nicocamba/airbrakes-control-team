%% ALTITUDE LINEAR MPC CONTROLLER. SUNSPEAR ROCKET. FILTERED SIGNAL
clc;
clear;
close all;

%% READING AND LOADING FLIGHT DATA
flight_structure = load("flight_data.mat");
flight_data = flight_structure.flight_data;

% At t=7.22s, thrust off
start_pos = find(flight_data.Time_s > 7.22);
start_pos = start_pos(1);

sat_on = find(flight_data.Time_s >= 7.8);
sat_on = sat_on(1); 

gravity = 9.77; 

S = 9.5e-3;     % Sref of the rocket

%% READING AND LOADING AB's DATA
Sab = 4e-3;
ab_drag_coeff = 0.22;

%% Filtrated data for density and temperature --> Only in control interval
air_dens_alt = [flight_data.Density_alt(start_pos : sat_on)' 
    flight_data.Z_m(start_pos : sat_on)']';

CD_speed = [flight_data.Total_speed_m_s(start_pos : sat_on)' 
    flight_data.CD_flight(start_pos : sat_on)']';
%% Simulation and timestep
h=0.01;    % Sample time       

time = 0:h:(flight_data.Time_s(sat_on) - flight_data.Time_s(start_pos));
steps = length(time);

%% Trajectory used as nominal --> u = 0
x0ref = [2253.934; 595];
input = 0;
steps_num = length(time)+100;

[altitude_ref, velocity_ref] = compute_ConstInput_apog(x0ref, input, steps_num, air_dens_alt,...
    CD_speed, ab_drag_coeff, h);

[z_apog_nom, z_apog_ind] = max(altitude_ref);

%% INITIAL VARIABLES
% State vector and control value
X0 = x0ref + [5;2];    % +10 m deviation and +2m/s in speed
U0 = 0;

% Initialization
x = X0;

position_sol = [X0(1) zeros(1,steps)];
velocity_sol = [X0(2) zeros(1,steps)];

tiempo_total = 0;

tiempo_iter_i = zeros(steps, 1);


%% Initialization of the kalman filter
% We initialize the estimated response (adding some error)
h_0 = x(1) + 15;
v_0 = x(2) + 3;

% Sensors characteristics
%IMU
BW = 10;
dens_noise_IMU = 230; % Technical specifications: micro g/sqrt(Hz)
var_IMU = (dens_noise_IMU*1e-6*gravity*sqrt(BW))^2;

%BIAS
var_bias = (50e-03*gravity/3)^2;
var_bias0 = 1e-4;
bias_var=0.1;

y_estim = [[h_0; v_0; 0] zeros(3,steps)]; % Initial deviation for filter

kalman_filter = Kalman_develop;
kalman_filter.var_IMU = var_IMU;
kalman_filter.var_bias0 = var_bias0;
kalman_filter.x_0 = [h_0; v_0; 0];
kalman_filter.P_0 = [15^2, 0, 0; 
                    0, 5^2, 0;
                    0, 0, var_bias];      % P0 from Thrust off!
kalman_filter.Ts = h;

kalman_filter.init_filter;

%% Linealization parameters
mass = 16.4813;
Cd_ab = 0.23;

% Noise generator for perturbations
processNoiseStd = 0.2;
total_steps = steps;
noise_seq = processNoiseStd * randn(total_steps, 1);
air_dens = 78181; 

%% Development
for s=1:total_steps
    xprev = x;
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

   velocity_sol(s+1) = x(2);
   position_sol(s+1) = x(1);

    %% State estimation
    bias_var = bias_var + randn*sqrt(var_bias0);   % random walk (o bias constante)
    u_Kal = (x(2) - xprev(2))/h + randn * sqrt(var_IMU) + bias_var; 

    %u_Kal = thrust_curve(s)/mass_T(s) - gravity - 0.5 * air_dens_r * x(2)^2 * S * drag_coeff_r  / mass_T(s) + ...
    %    noise + randn * sqrt(var_IMU);
    kalman_filter.prediction_step(u_Kal);


    % As there's no measurement, no gain

    kalman_filter.x_bar = kalman_filter.x_hat;
    

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
sgtitle('Kalman Filter Estimation Error during Supersonic Burnout Phase, T_s=0.01s', ...
        'FontSize', 14, 'FontWeight','bold')

