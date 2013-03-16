function [labels, primal_energy, dual_energy, time, oracle_calls] = bfgsDual(unary, vertC, horC)
	addpath('hanso/');

	[K, N, M] = size(unary);
	pars.nvar = K * N * M;
	wrapper = gridDualWrapper(unary, vertC, horC);
	function [value, derivative] = minus_dual(lambda, pars)
		[value, derivative] = wrapper.dual(reshape(lambda, K, N * M));
		value = -value;
		derivative = -derivative(:);
	end
	pars.fgname = @minus_dual;
	options.prtlevel = 2;
	options.maxit = 30;
	options.x0 = zeros(K * N * M, 1);
	[lambda, ~] = hanso(pars, options);
	[~, primal_energy, dual_energy, time] = wrapper.getState();
	[~, ~, ~, labels] = gridDual(reshape(lambda, K, N * M), unary, vertC, horC);
	oracle_calls = 1:length(dual_energy);
end
