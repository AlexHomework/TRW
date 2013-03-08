function [dual_energy, grad, upper_energy, labels_first, labels_second] = gridDual(lambda, unary, vertC, horC)
	% Minimize dual energy (grid potts model)
	% lambda is dual variabels matrix for horizontal chains
	% -lambda -- for vertical chains
	% 
	% grad -- is a projected subgradient
	% 

	[K, N] = size(lambda);

	dual_unary = unary / 2;
	[localEnergy, wholeEnergy, labels_first] = horizontalChains(lambda, unary, dual_unary, vertC, horC);
	dual_energy = localEnergy;
	upper_energy = wholeEnergy;
	[localEnergy, wholeEnergy, labels_second] = verticalChains(-1 * lambda, unary, dual_unary, vertC, horC);
	dual_energy = dual_energy + localEnergy;
	upper_energy = min(wholeEnergy, upper_energy);

	grad = sparse(K, N);
	for p = 1:K
		grad(p, :) = reshape(((labels_first == p) - (labels_second == p)), 1, N);
	end
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
