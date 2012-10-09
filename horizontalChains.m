
% function [localEnergy, wholeEnergy, labels] = horizontalChains(unary, dual_unary, vertC, horC, metric)
% 	[K, N, M] = size(unary);
% 	% lambda_unary = reshape(lambda_unary, K, N, M);
% 	labels = zeros(N, M);
% 	localEnergy = 0;
% 	for chain_i = 1:N
% 		chain_unary = reshape(dual_unary(:, chain_i, :), K, M);
% 		[sub_en, labels(chain_i, :)] = minimize_chain(chain_unary, horC(chain_i, :), metric);
% 		localEnergy = localEnergy + sub_en;
% 	end
% 	wholeEnergy = gridEnergy(unary, vertC, horC, metric, labels);
% 	labels = reshape(labels, N * M, 1);
% end
function [localEnergy, wholeEnergy, labels] = horizontalChains(unary, dual_unary, vertC, horC, metric)
	[K, N, M] = size(unary);
	labels = zeros(N, M);
	localEnergy = 0;

	metric = reshape(metric, K, 1, K);
	rep_K = ones(1, K); % we will use this matrix as index inspite of using repmat
	rep_N = ones(1, N);

	% message_forward(:, i) it's message from i to i + 1
	message_forward = zeros(K, N, M - 1);
	message_forward_y = zeros(K, N, M - 1);

	[message_forward(:, :, 1), message_forward_y(:, :, 1)] = ...
					reshape(min(dual_unary(:, :, rep_K) + ...
					bsxfun(@times, reshape(horC(:, 1), 1, N), ...
								   metric(:, rep_N, :)), [], 1), N, K);
	for t = 2:M-1
		el = dual_unary(:, :, t) + message_forward(:, :, t - 1);
		[message_forward(:, t), message_forward_y(:, t)] = reshape(min(el(:, :, rep_K) + ...
					bsxfun(@times, reshape(horC(:, t), 1, N), ...
								   metric(:, rep_N, :)), [], 1), N, K)';
	end

	y = zeros(N, M);
	[energy_arr, y(:, M)] = min(message_forward(:, :, M - 1) + dual_unary(:, :, M));
	energy = sum(energy_arr);
	for t = M-1:-1:1
		y(:, t) = message_forward_y(y(:, t + 1), :, t);
	end
end