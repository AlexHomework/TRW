function [dataCost, neighbors, metric] = toClassicFormat(unary, vertC, horC)
	[K, N, M] = size(unary);
	dataCost = reshape(unary, K, N * M);

	neighbors = sparse(N * M, N * M);
	% Add vertical edges
	for i = 1:M
		first = 1 + (i - 1) * N; % first element in the i-th column
		last = first + N - 1;
		neighbors(first : last-1, first+1 : last) = diag(vertC(:, i));
	end
	% Add horizontal edges
	for i = 1:M-1
		first = 1 + (i - 1) * N; % first element in the i-th column
		last = first + N - 1;
		neighbors(first : last, first+N : last+N) = diag(horC(:, i));
	end

	metric = ones(K, K) - eye(K);
end
