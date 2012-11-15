function [labels, energy, lowerBound, time] = TRW_S(unary, vertC, horC)
	addpath('TRW_S/');

	[K, N, M] = size(unary);
	[dataCost, neighbors, metric] = toClassicFormat(unary, vertC, horC);
	t = cputime;
	[labels, energy, lowerBound] = mrfMinimizeMex(dataCost, neighbors, metric);
	time = cputime - t;
	labels = reshape(labels, N, M);
end
