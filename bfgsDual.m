function [labels, energy, lowerBound, time] = bfgsDual(unary, vertC, horC)
	addpath('hanso/');

	[K, N, M] = size(unary);
	pars.nvar = K * N * M;
	function [value, derivative] = minus_dual(lambda, pars)
		[value, derivative] = gridDual(reshape(lambda, K, N * M), unary, vertC, horC);
		value = -value;
		derivative = -derivative(:);
	end
	pars.fgname = @minus_dual;
	options.prtlevel = 2;
	options.cpumax = 90;
	t = cputime;
	[lambda, ~] = hanso(pars, options);
	time = cputime - t;
	[lowerBound, ~, energy, labels] = gridDual(reshape(lambda, K, N * M), unary, vertC, horC);
end
