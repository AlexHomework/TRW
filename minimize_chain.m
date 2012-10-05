function [energy, y] = minimize_chain(unary, C, metric)
% Binary potential between i and i+1 is C(i) * metric(i, i+1)

% Length of the chain
N = size(unary, 1);
% Possible y values: from 1 to K
K = size(metric, 1);

idx_matrix = ones(1, K); % we will use this matrix as index inspite of using repmat

% message_forward(i) it's message from i to i + 1
message_forward = zeros(N - 1, K);
message_forward_y = zeros(K, N - 1);

[values, message_forward_y(:, 1)] = min(unary(idx_matrix, :) + C(1) * metric, [], 2);
message_forward(1, :) = values';
for t = 2:N-1
	el = unary(t, :) + message_forward(t - 1, :);
	[values, message_forward_y(:, t)] = min(el(idx_matrix, :) + C(t) * metric, [], 2);
	message_forward(t, :) = values';
end

y = zeros(1, N);
[energy, y(N)] = min(message_forward(N - 1, :) + unary(N, :));
for t = N-1:-1:1
	y(t) = message_forward_y(y(t + 1), t);
end
end