function [best_h, dual_energy] = maxBacktracking(func, func_0, grad_times_direction, h_init, f_init)
	% grad_times_direction is just dot product of
	% optimization direction and the gradient in func(0)
	% 
	
	max_iter = 15;
	fact = 0.5;
	rho = 0.1;
	h = h_init;
	it = 1;
	dual_energy = f_init;

	while (it <= max_iter) & (dual_energy < (func_0 + rho * h * grad_times_direction))
		h = h * fact;
		it = it + 1;
		dual_energy = func(h);
	end
	best_h = h;
	if it > max_iter
		disp('Reached maximum number of iterations!');
	else
		disp('Approximate optimization succeeded.');
	end
end