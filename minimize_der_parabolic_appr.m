function x_min = minimize_der_parabolic_appr(x_1, f_d_1, x_2, f_d_2)
	% Minimum point of parabolic approximation given derivatives in two points
	x_min = (f_d_2 * x_1 - f_d_1 * x_2) / (f_d_2 - f_d_1);
end