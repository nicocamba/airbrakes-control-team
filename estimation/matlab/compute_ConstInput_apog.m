function [altitude, velocity] = compute_ConstInput_apog(x0, input, steps_num, air_dens_alt, drag_coef_speed, ...
    ab_drag_coeff, h)
    altitude = [x0(1) zeros(1,steps_num-1)];
    velocity = [x0(2) zeros(1,steps_num-1)];
    x=x0;
    for i=1:steps_num
        [air_dens, drag_coeff] = compute_current_parameters(air_dens_alt, drag_coef_speed,...
            altitude(i), velocity(i));
    
        t_span = h*[i i+1];
        [t, X] = ode45(@(t, X) compute_dynamics(X, input, air_dens, drag_coeff, ab_drag_coeff), t_span, x);
        x =  X(end,:)';
        altitude(i+1) = X(end,1);
        velocity(i+1) = X(end,2);
    
    end
    
end