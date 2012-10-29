function [labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC)
[K, N, M] = size(unary);
dual_unary = unary / 2;

f1 = @(lambda) horizontalChains(lambda, unary, dual_unary, vertC, horC);
f2 = @(lambda) verticalChains(lambda, unary, dual_unary, vertC, horC);
[labels, energy, lowerBound, time] = dualDecomposition(K, N * M, f1, f2);
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
	labels = zeros(N, M);
	localEnergy = 0;
	for chain_i = 1:N
		chain_unary = reshape(dual_unary(:, chain_i, :), K, M) + reshape(lambda_unary(:, chain_i, :), K, M);
		[sub_en, labels(chain_i, :)] = minimize_chain(chain_unary, horC(chain_i, :));
		localEnergy = localEnergy + sub_en;
	end
	wholeEnergy = gridEnergy(unary, vertC, horC, labels);
	labels = reshape(labels, N * M, 1);
end
function [localEnergy, wholeEnergy, labels] = verticalChains(lambda_unary, unary, dual_unary, vertC, horC)
	[K, N, M] = size(unary);
	lambda_unary = reshape(lambda_unary, K, N, M);
	labels = zeros(N, M);
	localEnergy = 0;
	for chain_i = 1:M
		chain_unary = reshape(dual_unary(:, :, chain_i), K, N) + reshape(lambda_unary(:, :, chain_i), K, N);
		[sub_en, labels(:, chain_i)] = minimize_chain(chain_unary, vertC(:, chain_i));
		localEnergy = localEnergy + sub_en;
	end
	wholeEnergy = gridEnergy(unary, vertC, horC, labels);
	labels = reshape(labels, N * M, 1);
end
