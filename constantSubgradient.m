function [first_diff, context, alpha_n] = ...
					constantSubgradient(labels_first, labels_second, ...
					lowerBound, best_dual_energy, K, N, iteration, context)
	alpha_n = context.step;

	first_diff = zeros(K, N);
	% Lambda projected subgradient maximization
	for p = 1:K
		first_diff(p, :) = alpha_n * reshape(((labels_first == p) - (labels_second == p)), 1, N);
	end
end