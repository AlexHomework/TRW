function [labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, metric)
[N, M, K] = size(unary);
hor_y = zeros(N, M);
lambda_unary_hor = zeros(N, M, K);
ver_y = zeros(N, M);
lambda_unary_ver = zeros(N, M, K);
dual_unary = unary / 2;

f1 = @(lambda) horizontalChains(lambda, unary, dual_unary, vertC, horC, metric);
f2 = @(lambda) verticalChains(lambda, unary, dual_unary, vertC, horC, metric);
[labels, energy, lowerBound, time] = dualDecomposition(N * M, K, f1, f2);
labels = reshape(labels, N, M);

showImage(labels);

figure;
plot(energy);
hold on;
plot(lowerBound, 'r');
figure;
plot(time);

end


function [localEnergy, wholeEnergy, labels] = horizontalChains(lambda_unary, unary, dual_unary, vertC, horC, metric)
	[N, M, K] = size(unary);
	lambda_unary = reshape(lambda_unary, N, M, K);
	labels = zeros(N, M);
	localEnergy = 0;
	for chain_i = 1:N
		chain_unary = reshape(dual_unary(chain_i, :, :), M, K) + reshape(lambda_unary(chain_i, :, :), M, K);
		[sub_en, labels(chain_i, :)] = minimize_chain(chain_unary, horC(chain_i, :), metric);
		localEnergy = localEnergy + sub_en;
	end
	wholeEnergy = gridEnergy(unary, vertC, horC, metric, labels);
	labels = reshape(labels, N * M, 1);
end
function [localEnergy, wholeEnergy, labels] = verticalChains(lambda_unary, unary, dual_unary, vertC, horC, metric)
	[N, M, K] = size(unary);
	lambda_unary = reshape(lambda_unary, N, M, K);
	labels = zeros(N, M);
	localEnergy = 0;
	for chain_i = 1:M
		chain_unary = reshape(dual_unary(:, chain_i, :), N, K) + reshape(lambda_unary(:, chain_i, :), N, K);
		[sub_en, labels(:, chain_i)] = minimize_chain(chain_unary, vertC(:, chain_i), metric);
		localEnergy = localEnergy + sub_en;
	end
	wholeEnergy = gridEnergy(unary, vertC, horC, metric, labels);
	labels = reshape(labels, N * M, 1);
end