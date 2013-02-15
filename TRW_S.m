function [labels, energy, lowerBound, time] = TRW_S(unary, vertC, horC)
	addpath('TRW_S/');

	[K, N, M] = size(unary);
	[dataCost, neighbors, ~] = toClassicFormat(unary, vertC, horC);
	[labels, energy, lowerBound, time] = mrfMinimizeMex(dataCost, neighbors, []);
	% [labels, energy, lowerBound, time] = mrfMinimizeMex(dataCost, neighbors, [], struct('funcEps', -1, 'maxIter', 450));
	labels = reshape(labels, N, M);
end
