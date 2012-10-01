function [unary, vertC, horC, metric] = potentials(image1, image2)
% input: two YUV stereo images
[N, M, ~] = size(image1);
% Let's for begining fix disparities range to [-20, 20]
K = 41;
L = 20;
alpha = 1;
s = 2;

unary = zeros(N, M, K);
image1_sq = sum(image1.^2, 3);
image2_sq = sum(image2.^2, 3);
for p = 1:K
	dspr = disparity(p, K);
	if dspr == 0
		unary(:, :, p) = sqrt(image1_sq(:, :) + image2_sq(:, :));
	elseif dspr < 0
		dspr = -dspr;
		unary(:, dspr+1 : M, p) = sqrt(image1_sq(:, dspr+1 : M) + image2_sq(:, 1 : M-dspr));
		unary(:, 1 : dspr, p) = +inf * ones(N, dspr, 1);
	else
		unary(:, 1 : M-dspr, p) = sqrt(image1_sq(:, 1 : M-dspr) + image2_sq(:, dspr+1 : M));
		unary(:, M-dspr+1 : M, p) = +inf * ones(N, dspr, 1);
	end
end


vertC = alpha * exp(-abs(image1(1 : N-1, :) - image1(2 : N, :)) / s);
horC = alpha * exp(-abs(image1(:, 1 : M-1) - image1(:, 2 : M)) / s);


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