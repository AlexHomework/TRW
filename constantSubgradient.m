function [lambda_first_diff, lambda_second_diff, context] = ...
					constantSubgradient(labels_first, labels_second, ...
					lowerBound, best_dual_energy, K, N, iteration, context)
	alpha_n = 1;

	% Lambda projected subgradient maximization
	for p = 1:K
		lambda_first_diff(p, :) = alpha_n * reshape(((labels_first == p) - (labels_second == p)), 1, N);
		lambda_second_diff(p, :) = alpha_n * reshape(((labels_second == p) - (labels_first == p)), 1, N);
	end
end