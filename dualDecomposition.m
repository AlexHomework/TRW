function [labels, energy, lowerBound, time, step, dual_calls] = dualDecomposition(K, N, u_dual_func, ...
																dual_step, init_context, varargin)
	% Primal problem:
	% min f1'(x1, \theta1) + min f2'(x2, \theta2)
	% subject to x1 = x2
	% 
	% f1(\lambda) = min over x1 (f1'(x1, \theta1 + \lambda))
	% f2(\lambda) = min over x2 (f2'(x2, \theta2 - \lambda))
	% Dual problem:
	% find \lambda which maximize (f1(\lambda) + f2(\lambda))
	% 
	% \theta is a matrix:
	% [K, N] = size(\theta)
	% 
	% [dual_energy, grad, upper_energy, labels_first, labels_second] = dual_func(\lambda)
	% where:
	% dual_energy = f1(\lambda) + f2(\lambda)
	% grad is projected subgradient
	% upper_energy is current primal energy estimate
	% labels_first  = argmin over x1 (f1'(x1, \theta1 + \lambda))
	% labels_second = argmin over x2 (f2'(x2, \theta2 - \lambda))
	%
	% Optional parameter lambda is for lambda_first initialization
	% 

	[iterations_count, lambda_first] = process_options(varargin, 'iter', 600, 'lambda', zeros(K, N));
	if prod(size(lambda_first)) == 0
		lambda_first = zeros(K, N);
	end

	curr_dual_calls = 0;
	function [dual_energy, grad, upper_energy, labels_first, labels_second] = dual_func(lambda)
		curr_dual_calls = curr_dual_calls + 1;
		[dual_energy, grad, upper_energy, labels_first, labels_second] = u_dual_func(lambda);
	end

	lowerBound = [];
	energy = [];
	best_dual_energy = 0;
	time = [];
	step = [];
	dual_calls = [];
	context = init_context;
	t = cputime;
	for iteration = 1:iterations_count
		% We will skip this part in time counting because it was already
		% computed on the previous iteration in dual_step computation
		% but it's really hard to store it (using neat code), so we recompute it
		skip_time_start = cputime;
		skip_oracle_calls_start = curr_dual_calls;
		[dual_energy, grad, upper_energy, labels_first, labels_second] = dual_func(lambda_first);
		curr_dual_calls = skip_oracle_calls_start;
		t = t + (cputime - skip_time_start);

		best_dual_energy = max(best_dual_energy, dual_energy);
		lowerBound = [lowerBound; dual_energy];
		energy = [energy; upper_energy];


		[context, curr_step, curr_dual_energy] = dual_step(@(step) dual_func(lambda_first + step * grad), ...
													grad(:), lowerBound, iteration, context);

		dual_calls = [dual_calls; curr_dual_calls];
		step = [step; curr_step];

		lambda_first = lambda_first + curr_step * grad;
		time = [time; cputime - t];
	end

	labels = labels_first;
end
