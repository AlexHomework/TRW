function [localEnergy, labels] = minimizeHorChains(unary, horC)
	[K, N, M] = size(unary);

	rep_K = ones(1, K); % we will use this matrix as index inspite of using repmat
	rep_N = ones(1, N);

	% message_forward(:, i) it's message from i to i + 1
	message_forward = zeros(K, N, M - 1);
	message_forward_y = zeros(K, N, M - 1);

	values = unary(:, :, 1);
	[min_values, min_idx] = min(values, [], 1);
	min_values = min_values + horC(:, 1)';
	min_values = min_values(rep_K, :);
	min_idx = min_idx(rep_K, :);
	greater = values > min_values;
	values(greater) = min_values(greater);
	message_forward(:, :, 1) = values;
	indexes = (1:K)';
	indexes = indexes(:, rep_N);
	indexes(greater) = min_idx(greater);
	message_forward_y(:, :, 1) = indexes;
	for t = 2:M-1
		values = unary(:, :, t) + message_forward(:, :, t - 1);
		[min_values, min_idx] = min(values, [], 1);
		min_values = min_values + horC(:, t)';
		min_values = min_values(rep_K, :);
		min_idx = min_idx(rep_K, :);
		greater = values > min_values;
		values(greater) = min_values(greater);
		message_forward(:, :, t) = values;
		indexes = (1:K)';
		indexes = indexes(:, rep_N);
		indexes(greater) = min_idx(greater);
		message_forward_y(:, :, t) = indexes;
	end

	labels = zeros(N, M);
	[energy_arr, labels(:, M)] = min(message_forward(:, :, M - 1) + unary(:, :, M));
	localEnergy = sum(energy_arr);
	for t = M-1:-1:1
		labels(:, t) = message_forward_y(labels(:, t + 1)' + (0:(N - 1))*K + (t - 1)*K*N);
	end
end