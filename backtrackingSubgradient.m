function [context, alpha_n, f_n] = backtrackingSubgradient(func, direction, grad, ...
														lower_bound, iteration, context)
	% 1-d step optimization using backtracking method
	if isfield(context, 'use_adaptive_init') & context.use_adaptive_init
		[context, init, f_init] = adaptiveSubgradient(func, direction, grad, ...
												lower_bound, iteration, context);
	else
		if iteration == 1
			init = 1;
		else
			init = context.last_step;
		end
		f_init = func(init);
	end
	[alpha_n, f_n] = maxBacktracking(func, lower_bound(end), grad' * direction, init, f_init);
	context.last_step = alpha_n;
end