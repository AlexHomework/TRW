function [context, alpha_n, f_n] = adaptiveSubgradient(func, grad, lower_bound, iteration, context)
	 % Adaptive projected subgradient step computation
	gamma0 = 1.5;
	gamma1 = 0.5;
	epsilon = @(n) 1 / n;
	init_delta_prev = 1000;

	if iteration == 1
		delta = init_delta_prev;
	else
		if lower_bound(iteration) > lower_bound(iteration - 1)
			delta = gamma0 * context.delta_prev;
		else
			delta = max(gamma1 * context.delta_prev, epsilon(iteration));
		end
	end
	context.delta_prev = delta;
	best_dual_energy = max(lower_bound);
	alpha_n = best_dual_energy + delta - lower_bound(iteration);
	alpha_n = alpha_n / sum(abs(grad));

	f_n = func(alpha_n);
end