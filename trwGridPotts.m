function [labels, energy, lowerBound, time, step, dual_calls] = trwGridPotts(unary, vertC, horC, ...
												dual_step, init_context, varargin)
	[K, N, M] = size(unary);
	dual_unary = unary / 2;
	[random_init] = process_options(varargin, 'random_init', 0);

	if random_init
		% Generate random lambda initialization
		lambda = randn(K, N * M) / 2;
	else
		lambda = [];
	end
	
	[labels, energy, lowerBound, time, step, dual_calls] = dualDecomposition(K, N * M, ...
						@(lambda) gridDual(lambda, unary, vertC, horC), dual_step, ...
												init_context, 'lambda', lambda);
	
	labels = reshape(labels, N, M);

end
