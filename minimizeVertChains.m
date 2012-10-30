function [localEnergy, labels] = minimizeVertChains(unary, vertC)
	[K, N, M] = size(unary);
	labels = zeros(N, M);
	localEnergy = 0;

	rep_K = ones(1, K); % we will use this matrix as index inspite of using repmat
	rep_M = ones(1, M);

	% message_forward(:, i) it's message from i to i + 1
	message_forward = zeros(K, N - 1, M);
	message_forward_y = zeros(K, N - 1, M);

	values = reshape(unary(:, 1, :), K, M);
	[min_values, min_idx] = min(values, [], 1);
	min_values = min_values + vertC(1, :);
	min_values = min_values(rep_K, :);
	min_idx = min_idx(rep_K, :);
	greater = values > min_values;
	values(greater) = min_values(greater);
	message_forward(:, 1, :) = values;
	indexes = (1:K)';
	indexes = indexes(:, rep_M);
	indexes(greater) = min_idx(greater);
	message_forward_y(:, 1, :) = indexes;
	for t = 2:N-1
		values = reshape(unary(:, t, :) + message_forward(:, t - 1, :), K, M);
		[min_values, min_idx] = min(values, [], 1);
		min_values = reshape(min_values, 1, M);
		min_idx = reshape(min_idx, 1, M);
		min_values = min_values + vertC(t, :);
		min_values = min_values(rep_K, :);
		min_idx = min_idx(rep_K, :);
		greater = values > min_values;
		values(greater) = min_values(greater);
		message_forward(:, t, :) = values;
		indexes = (1:K)';
		indexes = indexes(:, rep_M);
		indexes(greater) = min_idx(greater);
		message_forward_y(:, t, :) = indexes;
	end

	labels = zeros(N, M);
	[energy_arr, labels(N, :)] = min(message_forward(:, N - 1, :) + unary(:, N, :));
	localEnergy = sum(energy_arr);
	for t = N-1:-1:1
		labels(t, :) = message_forward_y(labels(t + 1, :) + (t - 1)*K + (0:(M - 1))*K*(N-1));
	end
end