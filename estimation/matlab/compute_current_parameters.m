function [air_dens, drag_coeff] = compute_current_parameters(air_dens_alt, drag_coef_speed,  current_alt, current_vel)

    [~, air_dens_pos] = min(abs(air_dens_alt(2, :) - current_alt));
    air_dens = air_dens_alt(1, air_dens_pos);

    [~, drag_coeff_pos] = min(abs(drag_coef_speed(1, :) - current_vel));
    drag_coeff = drag_coef_speed(2, drag_coeff_pos);


end
