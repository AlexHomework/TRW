function [energy, y] = minimize_chain(unary, C)
% Binary potential between i and i+1 is C(i) * [y(i) ~= y(i + 1)]

% N is the length of the chain
% Possible y values: from 1 to K
[K, N] = size(unary);

idx_matrix = ones(1, K); % we will use this matrix as index inspite of using repmat

% message_forward(:, i) it's message from i to i + 1
message_forward = zeros(K, N - 1);
message_forward_y = zeros(K, N - 1);


values = unary(:, 1);
[min_value, min_idx] = min(values);
min_value = min_value + C(1);
idx = (values > min_value);
values(idx) = min_value;
indexes = 1:K;
indexes(idx) = min_idx;
message_forward(:, 1) = values;
message_forward_y(:, 1) = indexes;
for t = 2:N-1
	values = unary(:, t) + message_forward(:, t - 1);
	[min_value, min_idx] = min(values);
	min_value = min_value + C(t);
	idx = (values > min_value);
	values(idx) = min_value;
	indexes = 1:K;
	indexes(idx) = min_idx;
	message_forward(:, t) = values;
	message_forward_y(:, t) = indexes;
end

y = zeros(1, N);
[energy, y(N)] = min(message_forward(:, N - 1) + unary(:, N));
for t = N-1:-1:1
	y(t) = message_forward_y(y(t + 1), t);
end
end
