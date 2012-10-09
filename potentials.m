function [unary, vertC, horC, metric] = potentials(image1, image2)
% input: two YUV stereo images
[N, M, ~] = size(image1);
% Let's for begining fix disparities range to [-20, 20]
K = 61;
L = 10;
alpha = 1;
s = 7;

unary = zeros(K, N, M);
image1_sq = sum(image1.^2, 3);
image2_sq = sum(image2.^2, 3);
for p = 1:K
	dspr = disparity(p, K);
	if dspr == 0
		unary(p, :, :) = sqrt(sum((image1(:, :, :) - image2(:, :, :)) .^ 2, 3));
	elseif dspr < 0
		dspr = -dspr;
		unary(p, :, dspr+1 : M) = sqrt(sum((image1(:, dspr+1 : M, :) - image2(:, 1 : M-dspr, :)) .^ 2, 3));
		unary(p, :, 1 : dspr) = 1000 * ones(N, dspr, 1);
	else
		unary(p, :, 1 : M-dspr) = sqrt(sum((image1(:, 1 : M-dspr, :) - image2(:, dspr+1 : M, :)) .^ 2, 3));
		unary(p, :, M-dspr+1 : M) = 1000 * ones(N, dspr, 1);
	end
end


vertC = alpha * exp(-abs(image1(1 : N-1, :, 1) - image1(2 : N, :, 1)) / s);
horC = alpha * exp(-abs(image1(:, 1 : M-1, 1) - image1(:, 2 : M, 1)) / s);


metric = L * ones(K, K);
for i = 1:K
	for j = i:min(K, i + L - 1)
		metric(i, j) = j - i;
	end
	for j = max(1, i - L + 1):i
		metric(i, j) = i - j;
	end
end
end


function disparity = disparity(p, K)
	disparity = p - ceil(K / 2);
end
