function [best_h, dual_energy] = maxBacktracking(func, func_0, grad_sq, h_init)
	max_iter = 15;
	fact = 0.5;
	rho = 0.1;
	h = h_init;
	it = 1;
	dual_energy = func(h);

	while (it <= max_iter) & (dual_energy < (func_0 + rho * h * grad_sq))
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