function thrust_curve = thrust_interpolation(thrust_data, h)

t = thrust_data(:,1);
F = thrust_data(:,2);

t_new = (0:h:t(end))';      % vector columna

F_new = interp1(t, F, t_new, 'linear');

thrust_curve = [t_new, F_new];

end