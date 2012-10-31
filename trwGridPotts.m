function [labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC)
[K, N, M] = size(unary);
dual_unary = unary / 2;

f1 = @(lambda) horizontalChains(lambda, unary, dual_unary, vertC, horC);
f2 = @(lambda) verticalChains(lambda, unary, dual_unary, vertC, horC);
[labels, energy, lowerBound, time] = dualDecomposition(K, N * M, f1, f2, @adaptiveSubgradient);
labels = reshape(labels, N, M);

showImage(labels);

figure;
plot(energy);
hold on;
plot(lowerBound, 'r');
figure;
plot(time);

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
