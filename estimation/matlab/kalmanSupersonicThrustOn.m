%% ALTITUDE LINEAR MPC CONTROLLER. SUNSPEAR ROCKET. FILTERED SIGNAL
clc;
clear;
close all;

%% READING AND LOADING FLIGHT DATA
flight_structure = load("flight_data.mat");
flight_data = flight_structure.flight_data;

thrust_structure = load("thrust_data.mat");
thrust_data = thrust_structure.thrust_data;

% Supersonic entry
start_pos = find(flight_data.Time_s >=  2.67);
start_pos = start_pos(1);

thrustOff = find(flight_data.Time_s >=7.22);
thrustOff = thrustOff(1); 

gravity = 9.77; 

S = 9.5e-3;     % Sref of the rocket

%% Filtrated data for density and temperature --> Only in control interval
air_dens_alt = [flight_data.Density_alt(start_pos : thrustOff)' 
    flight_data.Z_m(start_pos : thrustOff)']';

CD_speed = [flight_data.Total_speed_m_s(start_pos : thrustOff)' 
    flight_data.CD_flight(start_pos : thrustOff)']';


%% Simulation and timestep
h=0.01;    % Sample time       

time = 0:h:(flight_data.Time_s(thrustOff) - flight_data.Time_s(start_pos));
steps = length(time);

%% Noise generator for perturbations
processNoiseStd = 0.4*sqrt(h/0.1);
total_steps = steps;
noise_seq = processNoiseStd * randn(total_steps, 1);


%% Thrust curve interpolation in subsonic time range
thrust_curve = thrust_interpolation(thrust_data, h);
thrust_starting_index = find(thrust_curve(:,1) >=  2.67);
thrust_starting_index = thrust_starting_index(1);
thrust_curve = thrust_curve(thrust_starting_index:thrust_starting_index+steps,2);


mass_T =[22.715; zeros(size(thrust_curve,1)-1,1)];
Isp = 170.7; g = 9.77;

for i=1:(size(thrust_curve,1)-1)
    mass_dot = thrust_curve(i)/(Isp*g);  
    mass_T(i+1) = mass_T(i) - mass_dot*h;  
end

%% INITIAL VARIABLES
% State vector and control value
X0 = [1282.5173;320];   

% Initialization
x = X0;
position_sol = [X0(1) zeros(1,steps)];
velocity_sol = [X0(2) zeros(1,steps)];

tiempo_total = 0;

tiempo_iter_i = zeros(steps, 1);


%% Initialization of the kalman filter
% We initialize the estimated response (adding some error)
h_0 = x(1)+5;
v_0 = x(2)+2;

% Sensors characteristics

%IMU
BW = 10;
dens_noise_IMU = 230; % Technical specifications: micro g/sqrt(Hz)
var_IMU = (dens_noise_IMU*1e-6*gravity*sqrt(BW))^2;

%BIAS
var_bias = (50e-03*gravity/3)^2;
var_bias0 = 1e-4;
bias_var=0.05;

y_estim = [[h_0; v_0; 0] zeros(3,steps)]; % Initial deviation for filter

kalman_filter = Kalman_develop;
kalman_filter.var_model = processNoiseStd^2;
kalman_filter.var_IMU = var_IMU;
kalman_filter.var_bias0 = var_bias0;
kalman_filter.x_0 = [h_0; v_0; 0];
kalman_filter.P_0 = [5^2, 0, 0; 
                    0, 2^2, 0;
                    0, 0, var_bias];      % P0 from last stage!!!
kalman_filter.Ts = h;
kalman_filter.init_filter;

%% Development
for s=1:total_steps
    if mass_T(s)<=16.4813
        mass_T(s) = 16.4813;
        thrust_curve(s) = 0;
    end
    xprev = x;
    %% Determination CD and rho --> Using y of estimation!
    [air_dens, drag_coeff] = compute_current_parameters(air_dens_alt, CD_speed,...
        y_estim(1,s), y_estim(2,s));
   
   %% Input computed --> Integration of the system
    % Determination CD and rho for real system
   [air_dens_r, drag_coeff_r] = compute_current_parameters(air_dens_alt, CD_speed,...
      x(1), x(2));

   t_span = h*[s s+1];

   noise = noise_seq(s);
   [t, Xode] = ode45(@(t, X) compute_dynamics_thrust(X, air_dens_r, drag_coeff_r,...
       thrust_curve(s), mass_T(s), noise), t_span, x);
   
   x = Xode(end,:)';

   %% Save real simulation
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

%% Plots


figure;
plot(velocity_sol, position_sol, 'r-', 'LineWidth', 1.5); hold on;
xlabel('v');
ylabel('z');
title(" SUPERSONIC THRUST ON STAGE ");
legend('z(m)');
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
sgtitle('Kalman Filter Estimation Error during Supersonic Thrust Phase, T_s=0.01s', ...
        'FontSize', 14, 'FontWeight','bold')
