function [labels, energy, lowerBound, time] = TRW_S(unary, vertC, horC)
	addpath('TRW_S/');

	[K, N, M] = size(unary);
	[dataCost, neighbors, metric] = toClassicFormat(unary, vertC, horC);
	[labels, energy, lowerBound, time] = mrfMinimizeMex(dataCost, neighbors, metric, struct('funcEps', -1, 'maxIter', 500));
	labels = reshape(labels, N, M);
end
