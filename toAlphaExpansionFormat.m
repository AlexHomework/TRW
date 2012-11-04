function h = toAlphaExpansionFormat(unary, vertC, horC)
	% Return alpha expansion object
	addpath('gco-v3.0/matlab');

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

	h = GCO_Create(N * M, K);
	GCO_SetDataCost(h, dataCost);
	GCO_SetSmoothCost(h, metric);
	GCO_SetNeighbors(h, neighbors);
end
