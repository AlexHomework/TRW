function [labels, energy, lowerBound, time, step, dual_calls, iterations_info] = trwGridPotts(...
											unary, vertC, horC, dual_step, init_context, varargin)
	[K, N, M] = size(unary);
	dual_unary = unary / 2;
	[random_init, save_iterations] = process_options(varargin, 'random_init', 0, 'save_iterations', []);

	if random_init
		% Generate random lambda initialization
		lambda = randn(K, N * M) / 2;
	else
		lambda = [];
	end
	
	[labels, energy, lowerBound, time, step, dual_calls, iterations_info] = dualDecomposition(K, N * M, ...
						@(lambda) gridDual(lambda, unary, vertC, horC), dual_step, ...
						init_context, 'lambda', lambda, 'save_iterations', save_iterations);
	
	labels = reshape(labels, N, M);

end
