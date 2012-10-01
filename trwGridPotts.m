function [labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, metric)
[N, M, K] = size(unary);

hor_y = zeros(N, M);
lambda_unary_hor = zeros(N, M, K);
ver_y = zeros(N, M);
lambda_unary_ver = zeros(N, M, K);
lower_energy_arr = [];
upper_energy_arr = [];
for iteration = 1:10
	alpha_n = 0.5;

	% Y minimization
	% The lower energy estimate
	lower_energy = 0;
	% Horizontal chains
	for chain_i = 1:N
		chain_unary = reshape(unary(chain_i, :, :), M, K) + reshape(lambda_unary_hor(chain_i, :, :), M, K);
		[sub_en, hor_y(chain_i, :)] = minimize_chain(chain_unary, horC(chain_i, :), metric);
		lower_energy = lower_energy + sub_en;
	end

	% Vertical chains
	for chain_i = 1:M
		chain_unary = reshape(unary(:, chain_i, :), N, K) + reshape(lambda_unary_ver(:, chain_i, :), N, K);
		[sub_en, ver_y(:, chain_i)] = minimize_chain(chain_unary, vertC(:, chain_i), metric);
		lower_energy = lower_energy + sub_en;
	end

	upper_energy = min(gridEnergy(unary, vertC, horC, metric, hor_y), gridEnergy(unary, vertC, horC, metric, ver_y));

	lower_energy_arr = [lower_energy_arr, lower_energy];
	upper_energy_arr = [upper_energy_arr, upper_energy];

	% % Lambda subgradient maximisation
	for p = 1:K
		lambda_unary_hor(:, :, p) = lambda_unary_hor(:, :, p) + alpha_n * ((hor_y == p) - (ver_y == p));
		lambda_unary_ver(:, :, p) = lambda_unary_ver(:, :, p) + alpha_n * ((ver_y == p) - (hor_y == p));
	end
	
end

plot(upper_energy_arr);
hold on;
plot(lower_energy_arr, 'r');

end