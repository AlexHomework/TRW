function h = toAlphaExpansionFormat(unary, vertC, horC)
	% Return alpha expansion object
	addpath('alpha_expansion/matlab');

	[K, N, M] = size(unary);
	[dataCost, neighbors, metric] = toClassicFormat(unary, vertC, horC);

	h = GCO_Create(N * M, K);
	GCO_SetDataCost(h, dataCost);
	GCO_SetSmoothCost(h, metric);
	GCO_SetNeighbors(h, neighbors);
end
