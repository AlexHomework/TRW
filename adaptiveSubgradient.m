function [first_diff, context, alpha_n] = ...
						adaptiveSubgradient(labels_first, labels_second, lowerBound, ...
						best_dual_energy, K, N, iteration, context)
	 % Adaptive projected subgradient step computation
	gamma0 = 1.5;
	gamma1 = 0.5;
	epsilon = @(n) 1 / n;
	init_delta_prev = 1000;

	if iteration == 1
		delta = init_delta_prev;
	else
		if lowerBound(iteration) > lowerBound(iteration - 1)
			delta = gamma0 * context.delta_prev;
		else
			delta = max(gamma1 * context.delta_prev, epsilon(iteration));
		end
	end
	context.delta_prev = delta;
	alpha_n = best_dual_energy + delta - lowerBound(iteration);
	alpha_n = alpha_n / sum(sum(labels_first ~= labels_second));

	first_diff = zeros(K, N);
	% Lambda projected subgradient maximization
	for p = 1:K
		first_diff(p, :) = alpha_n * reshape(((labels_first == p) - (labels_second == p)), 1, N);
	end
end