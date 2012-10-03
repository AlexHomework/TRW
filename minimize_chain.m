function [energy, y] = minimize_chain(unary, C, metric)
% Binary potential between i and i+1 is C(i) * metric(i, i+1)

% Length of the chain
N = size(unary, 1);
% Possible y values: from 1 to K
K = size(metric, 1);

idx_matrix = ones(1, K); % we will use this matrix as index inspite of using repmat

% message_forward(i) it's message from i to i + 1
message_forward = zeros(N - 1, K);

message_forward(1, :) = min(unary(idx_matrix, :) + C(1) * metric, [], 2)';
for t = 2:N-1
	el = unary(t, :) + message_forward(t - 1, :);
	message_forward(t, :) = min(el(idx_matrix, :) + C(t) * metric, [], 2)';
end

% message_backward(i) it's message from i to i - 1
message_backward = zeros(N, K);

message_backward(N, :) = min(unary(N * idx_matrix, :) + C(N - 1) * metric, [], 2)';
for t = N-1:-1:2
	el = unary(t, :) + message_backward(t + 1, :);
	message_backward(t, :) = min(el(idx_matrix, :) + C(t - 1) * metric, [], 2)';
end

y = zeros(1, N);
[energy, y(1)] = min(message_backward(2, :) + unary(1, :));
[~, y(N)] = min(message_forward(N - 1, :) + unary(N, :));
for t = 2:N-1
	[~, y(t)] = min(message_forward(t - 1, :) + message_backward(t + 1, :) + unary(t, :));
end
end