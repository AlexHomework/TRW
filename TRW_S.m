function [labels, energy, lowerBound, time] = TRW_S(unary, vertC, horC, varargin)
	addpath('TRW_S/');
	[max_iter] = process_options(varargin, 'maxIter', 100);

	[K, N, M] = size(unary);
	[dataCost, neighbors, ~] = toClassicFormat(unary, vertC, horC);
	[labels, energy, lowerBound, time] = mrfMinimizeMex(dataCost, neighbors, [], struct('funcEps', -1, 'maxIter', max_iter));
	labels = reshape(labels, N, M);
end
