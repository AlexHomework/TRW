function [context, alpha_n, f_n] = constantSubgradient(func, grad, lower_bound, iteration, context)
	alpha_n = context.step;
	f_n = func(alpha_n);
end