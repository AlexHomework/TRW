function [context, alpha_n, f_n] = accurateSubgradient(func, grad, lower_bound, iteration, context)
	% Adaptive projected subgradient step computation
	% using exact (not approximate) optimization methods.
	
	[alpha_n, f_n] = fminbnd(@(x) -func(x), 0, 1.5);
	f_n = -f_n;
end