%% First Kalman model
classdef Kalman_develop < handle

    properties
        % Physical model variance
        var_model
        % Sensors variances
        var_IMU 
        var_baro 
        var_bias0
        var_speed

        % Sample time
        Ts

        % Current and predicted state estimations
        x_0 
        x_hat       % Future state prediction
        x_bar       % Current state estimation

        % Model matrices
        Ad
        Bd
        G
        H = [1, 0, 0];  % No changes

        % Covariance matrices of noise measurement
        Q % To be completed
        R

        % Kalman gain initialization
        K = zeros(3, 1);

        % Covariance matrix + initialization covariance matrix of state
        % estimation
        P_0
        P_hat
        P_bar
        
    end

    methods

        function set.Q(obj, Q)
            obj.Q = Q;
        end

        function set.R(obj, R)
            obj.R = R;
        end

        function set.var_model(obj, var_model)
            obj.var_model = var_model;
        end

        function set.var_speed(obj, var_speed)
            obj.var_speed = var_speed;
        end

        function set.var_IMU(obj, var_IMU)
            obj.var_IMU = var_IMU;
        end

        function set.var_bias0(obj, var_bias0)
            obj.var_bias0 = var_bias0;
        end
        
        function set.var_baro(obj, var_baro)
            obj.var_baro = var_baro;
        end

        function set.x_0(obj, x_0)
            obj.x_0 = x_0;
        end

        function set.P_0(obj, P_0)
            obj.P_0 = P_0;
        end

        function init_filter(obj)
            obj.x_bar = obj.x_0;
            obj.P_bar = obj.P_0;
            obj.Q = blkdiag((obj.var_IMU + obj.var_model), obj.var_bias0); % (Process noise, bias) 
            obj.Ad = [1, obj.Ts, -obj.Ts^2/2;
                0, 1, -obj.Ts;
                0, 0, 1];
            obj.Bd = [obj.Ts ^ 2/ 2;
                obj.Ts;
                0];
            obj.G = [obj.Ts ^ 2 / 2, 0;
                obj.Ts, 0;
                0, 1];

        end

        function prediction_step(obj, measured_u)
            obj.x_hat = obj.Ad * obj.x_bar + obj.Bd * measured_u;
            obj.P_hat = obj.Ad * obj.P_bar * obj.Ad' + obj.G * obj.Q * obj.G';
        end

        function update_step(obj, measured_y)
            obj.K = obj.P_hat * obj.H' / (obj.H * obj.P_hat * obj.H' + obj.R);
            obj.x_bar = obj.x_hat + obj.K * (measured_y - obj.H * obj.x_hat);
            obj.P_bar = (eye(3) - obj.K * obj.H) * obj.P_hat;
        end

    end

end