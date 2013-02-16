function [labels, energy, lowerBound, time, step, dual_calls, iterations_info] = ...
																dualDecomposition(K, N, u_dual_func, ...
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
	% Optional parameters are:
	%	lambda for lambda_first initialization
	% 	save_iterations is vector with iteration numbers,
	% 		about which dualDecomposition will return all possible
	% 		information (such as optimization direction)
	% 		in iterations_info vector

	[iterations_count, lambda_first, save_iterations] = process_options(varargin, 'iter', 400, ...
												'lambda', zeros(K, N), 'save_iterations', []);
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
	iterations_info = {};
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
		
		if (any(save_iterations == iteration))
			iterations_info{end + 1} = struct('iteration', iteration, 'step', curr_step, ...
											'oracle_calls', curr_dual_calls - dual_calls(end), ...
											'upper_energy', upper_energy, 'dual_energy', dual_energy, ...
											'direction', grad, 'lambda_first', lambda_first, ...
											'labels_first', labels_first, 'labels_second', labels_second);
		end

		dual_calls = [dual_calls; curr_dual_calls];
		step = [step; curr_step];

		lambda_first = lambda_first + curr_step * grad;
		time = [time; cputime - t];
	end

	labels = labels_first;
end
