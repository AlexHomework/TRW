function [context, alpha_n, f_n] = fletcherSubgradient(func, direction, grad, lower_bound, iteration, context)
	%  1-d step optimization using Fletcher's method

	[alpha_n, f_n] = maxFletcher(func, direction, lower_bound(end), grad, 'display', 1);
end