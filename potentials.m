function [unary, vertC, horC] = potentials(left, right)
% input: two YUV stereo images
[N, M, ~] = size(left);
K = 16;
alpha = 1;
s = 7;

unary = zeros(K, N, M);
for p = 1:K
	dspr = p - 1;
	if dspr > 0
		unary(p, :, dspr+1 : M) = sqrt(sum((left(:, dspr+1 : M, :) - right(:, 1 : M-dspr, :)) .^ 2, 3));
		% Some huge energies for impossible disparities
		unary(p, :, 1 : dspr) = 1000 * ones(N, dspr, 1);
	else
		unary(p, :, :) = sqrt(sum((left(:, :, :) - right(:, :, :)) .^ 2, 3));
	end
end


% vertC = alpha * exp(-abs(left(1 : N-1, :, 1) - left(2 : N, :, 1)) / s);
% horC = alpha * exp(-abs(left(:, 1 : M-1, 1) - left(:, 2 : M, 1)) / s);
vertC = 20 * ones(N - 1, M);
horC = 20 * ones(N, M - 1);
