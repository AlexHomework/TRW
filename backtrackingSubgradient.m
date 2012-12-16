function [context, alpha_n, f_n] = backtrackingSubgradient(func, grad, lower_bound, iteration, context)
	% Adaptive projected subgradient step computation using backtracking method
	if isfield(context, 'use_adaptive_init') & context.use_adaptive_init
		[context, init] = adaptiveSubgradient(func, grad, lower_bound, iteration, context);
	else
		if iteration == 1
			init = 1;
		else
			init = context.last_step;
		end
	end
	[alpha_n, f_n] = maxBacktracking(func, lower_bound(end), sum(abs(grad)), init);
	context.last_step = alpha_n;
end