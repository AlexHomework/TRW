function [labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, metric)
gamma0 = 1.5;
gamma1 = 0.5;
epsilon = @(n) 1 / n;
delta_prev = 1000;

[N, M, K] = size(unary);
hor_y = zeros(N, M);
lambda_unary_hor = zeros(N, M, K);
ver_y = zeros(N, M);
lambda_unary_ver = zeros(N, M, K);
dual_unary = unary / 2;
dual_energy_arr = [];
upper_energy_arr = [];
best_dual_energy = 0;
t = cputime;
for iteration = 1:10
	% Y minimization
	% The lower energy estimate
	dual_energy = 0;
	% Horizontal chains
	for chain_i = 1:N
		chain_unary = reshape(dual_unary(chain_i, :, :), M, K) + reshape(lambda_unary_hor(chain_i, :, :), M, K);
		[sub_en, hor_y(chain_i, :)] = minimize_chain(chain_unary, horC(chain_i, :), metric);
		dual_energy = dual_energy + sub_en;
	end

	figure;
	imshow(mat2gray(hor_y));

	% Vertical chains
	for chain_i = 1:M
		chain_unary = reshape(dual_unary(:, chain_i, :), N, K) + reshape(lambda_unary_ver(:, chain_i, :), N, K);
		[sub_en, ver_y(:, chain_i)] = minimize_chain(chain_unary, vertC(:, chain_i), metric);
		dual_energy = dual_energy + sub_en;
	end
	% figure;
	% imshow(mat2gray(ver_y));

	upper_energy = min(gridEnergy(unary, vertC, horC, metric, hor_y), gridEnergy(unary, vertC, horC, metric, ver_y));

	best_dual_energy = max(best_dual_energy, dual_energy);
	dual_energy_arr = [dual_energy_arr, dual_energy];
	upper_energy_arr = [upper_energy_arr, upper_energy];


	% Subgradient step computation
	if iteration == 1
		delta = delta_prev;
	else
		if dual_energy > dual_energy_arr(iteration - 1)
			delta = gamma0 * delta_prev;
		else
			delta = max(gamma1 * delta_prev, epsilon(iteration));
		end
	end
	delta_prev = delta
	alpha_n = best_dual_energy + delta - dual_energy
	best_dual_energy
	dual_energy
	alpha_n = alpha_n / sum(sum(hor_y ~= ver_y))
	% alpha_n = 0.5;


	% Lambda subgradient maximisation
	for p = 1:K
		lambda_unary_hor(:, :, p) = lambda_unary_hor(:, :, p) + alpha_n * ((hor_y == p) - (ver_y == p));
		lambda_unary_ver(:, :, p) = lambda_unary_ver(:, :, p) + alpha_n * ((ver_y == p) - (hor_y == p));
	end

	cputime - t
end

figure;
plot(upper_energy_arr);
hold on;
plot(dual_energy_arr, 'r');

end