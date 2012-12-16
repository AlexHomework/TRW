function [context, alpha_n, f_n] = fletcherSubgradient(func, grad, lower_bound, iteration, context)
	% Adaptive projected subgradient step computation using fletcher's method

	[alpha_n, f_n] = maxFletcher(func, grad, 'display', 1);
end