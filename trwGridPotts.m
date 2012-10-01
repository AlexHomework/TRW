function [labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, metric)

% chain_length = N * M;
% lambda_unary           = ones(chain_length, K);
% lambda_binary          = ones(chain_length - 1, K, K);

% chain_unary            = zeros(T, chain_length, K);
% % Horizontal snake, unary potentials
% for i = 1:N
% 	if (rem(i,2) == 0)
% 		% Even row
% 		chain_unary(1, (i - 1)*M + 1:i*M, :) = fliplr(unary(i, :));
% 	else
% 		% Odd row
% 		chain_unary(1, (i - 1)*M + 1:i*M, :) = unary(i, :);
% 	end
% end

% % Vertical snake, unary potentials
% for j = 1:M
% 	if (rem(j,2) == 0)
% 		% Even row
% 		chain_unary(1, (j - 1)*N + 1:j*N, :) = fliplr(unary(:, j));
% 	else
% 		% Odd row
% 		chain_unary(T, (j - 1)*N + 1:j*N, :) = unary(:, j);
% 	end
% end

% chain_binary           = zeros(T, chain_length - 1, K, K);
% % Horizontal snake, binary potentials
% for j = 1:M
% 	if (rem(j,2) == 0)
% 		% Even row
% 		chain_unary(1, (j - 1)*N + 1:j*N, :) = fliplr(unary(:, j));
% 	else
% 		% Odd row
% 		chain_unary(T, (j - 1)*N + 1:j*N, :) = unary(:, j);
% 	end
% end
% lambda_binary          = ones(chain_length - 1, K, K);
% lambda_binary(2, :, :) = -1 * lambda_binary;
% y = zeros(T, chain_length);

hor_y = zeros(N, M);

for iteration = 1:3
	alpha_n = 0.5;

	% Y minimization
	% The lower energy estimate
	% lower_energy = 0;
	% Horizontal chains
	for chain_i = 1:N
		[sub_en, hor_y(chain_i, :)] = minimize_chain(lambda_unary_hor(cain_i, :), horC(chain_i), metric);
		% lower_energy += sub_en;
	end

	% Vertical chains
	for chain_i = 1:M
		[sub_en, ver_y(:, chain_i)] = minimize_chain(lambda_unary_ver(:, chain_i), verC(chain_i), metric);
		% lower_energy += sub_en;
	end


	% Lambda subgradient maximisation
	lambda_unary_ver += alpha_n * (ver_y - hor_y);
	lambda_unary_hor += alpha_n * (hor_y - ver_y);
end
end