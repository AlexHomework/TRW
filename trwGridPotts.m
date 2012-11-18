function [labels, energy, lowerBound, time, step] = trwGridPotts(unary, vertC, horC, dualStep, ...
												init_context, varargin)
	[K, N, M] = size(unary);
	dual_unary = unary / 2;


	if length(varargin) > 0 && strcmp(varargin{1}, 'random_init')
		% Generate random lambda initialization
		lambda = randn(K, N * M) / 2;
	else
		lambda = [];
	end

	f1 = @(lambda) horizontalChains(lambda, unary, dual_unary, vertC, horC);
	f2 = @(lambda) verticalChains(lambda, unary, dual_unary, vertC, horC);
	[labels, energy, lowerBound, time, step] = dualDecomposition(K, N * M, f1, f2, dualStep, ...
													init_context, 'lambda', lambda);
	labels = reshape(labels, N, M);

end


function [localEnergy, wholeEnergy, labels] = horizontalChains(lambda_unary, unary, dual_unary, vertC, horC)
	[K, N, M] = size(unary);
	lambda_unary = reshape(lambda_unary, K, N, M);
	[localEnergy, labels] = minimizeHorChains(dual_unary + lambda_unary, horC);
	wholeEnergy = gridEnergy(unary, vertC, horC, labels);
end
function [localEnergy, wholeEnergy, labels] = verticalChains(lambda_unary, unary, dual_unary, vertC, horC)
	[K, N, M] = size(unary);
	lambda_unary = reshape(lambda_unary, K, N, M);
	[localEnergy, labels] = minimizeVertChains(dual_unary + lambda_unary, vertC);
	wholeEnergy = gridEnergy(unary, vertC, horC, labels);
end
