function [labels, energy, lowerBound, time, step, dual_calls] = subgradientDual(...
											unary, vertC, horC, dual_step, init_context, varargin)
	% Optimize dual energy via subgradient ascent.
	% Step 1-d optimization (step choosing) use dual_step function.
	% 

	[K, N, M] = size(unary);
	dual_unary = unary / 2;
	[random_init, max_iter, drawProfilePlot, profile_iterations, profileCallback] = ...
						process_options(varargin, 'randomInit', false, 'maxIter', 100, 'profile', false, ...
						'profileIters', [5, 10], 'profileCallback', @(fig, iter) 0);

	if drawProfilePlot
		disp(['Profile plots will be generated. Be careful, it', ...
			  ' can take a while and thus distort time measurements.']);
	end

	if random_init
		% Generate random lambda initialization
		lambda_first = randn(K, N * M) / 2;
	else
		lambda_first = zeros(K, N * M);
	end
	
	curr_dual_calls = 0;
	function [dual_energy, grad, upper_energy, labels_first, labels_second] = dual_func(lambda)
		curr_dual_calls = curr_dual_calls + 1;
		[dual_energy, grad, upper_energy, labels_first, labels_second] = gridDual(lambda, unary, vertC, horC);
	end

	lowerBound = zeros(max_iter, 1);
	energy = zeros(max_iter, 1);
	best_dual_energy = 0;
	time = zeros(max_iter, 1);
	step = zeros(max_iter, 1);
	dual_calls = zeros(max_iter, 1);
	context = init_context;
	t = cputime;
	for iteration = 1:max_iter
		% We will skip this part in time counting because it was already
		% computed on the previous iteration in dual_step computation
		% but it's really hard to store it (using neat code), so we recompute it
		skip_time_start = cputime;
		skip_oracle_calls_start = curr_dual_calls;
		[dual_energy, grad, upper_energy, labels_first, labels_second] = dual_func(lambda_first);
		curr_dual_calls = skip_oracle_calls_start;
		t = t + (cputime - skip_time_start);

		best_dual_energy = max(best_dual_energy, dual_energy);
		lowerBound(iteration) = dual_energy;
		energy(iteration) = upper_energy;

		direction = grad; % Compute optimization direction

		[context, curr_step, curr_dual_energy] = dual_step(@(step) dual_func(lambda_first + step * direction), ...
													direction(:), grad(:), lowerBound, iteration, context);

		if drawProfilePlot && any(iteration == profile_iterations)
			title = ['Profile subgradient (step size = ', int2str(curr_step)];
			title = [title, ', iter = ', int2str(iteration), ')'];
			profFig = profilePlot(@dual_func, lambda_first, lambda_first + curr_step * grad, 'title', title);
			profileCallback(profFig, iteration);
		end
		
		dual_calls(iteration) = curr_dual_calls;
		step(iteration) = curr_step;

		lambda_first = lambda_first + curr_step * grad;
		time(iteration) = cputime - t;

	end

	labels = reshape(labels_first, N, M);
end
