function [labels, energy, lowerBound, time, oracle_calls] = aggregateBundleDual(unary, vertC, horC, varargin)
	% Optimize dual energy via aggregate bundle method
	% 

	[K, N, M] = size(unary);
	[m, iterations_count, drawProfilePlot, profile_iterations, profileCallback] = ...
											process_options(varargin, ...
											'm', 0.1, 'maxIter', 40, 'profile', false, ...
											'profileIters', [5, 10], 'profileCallback', @(fig, iter) 0);


	precission = 1e-3;
	lambda_first = zeros(K * N * M, 1);

	if drawProfilePlot
		disp(['Profile plots will be generated. Be careful, it', ...
			  ' can take a while and thus distort time measurements.']);
	end

	wrapper = gridDualWrapper(unary, vertC, horC);
	function [dual_energy, grad, upper_energy, labels_first, labels_second] = dual(lambda)
		[dual_energy, grad, upper_energy, labels_first, labels_second] = ...
								wrapper.dual(reshape(lambda, K, N * M));
		% We want to maximize function using minimization algorithm
		grad = -grad(:);
		dual_energy = -dual_energy;
	end


	curr_x = lambda_first; % x_k
	prev_x = curr_x;
	[curr_x_dual_energy, curr_grad, upper_energy, labels_first] = dual(curr_x);
	prev_x_dual_energy = curr_x_dual_energy; % f(x_{k-1}) = f(x_k) for k = 1
	curr_aggregate = curr_grad;
	curr_epsilon = 0; % \epsilon_{k-1}
	curr_alpha = 0; % \alpha_k
	for iteration = 1:iterations_count
		curr_mu = fundDirection(curr_aggregate, curr_grad, curr_epsilon, curr_alpha);

		curr_aggregate = (1 - curr_mu) * curr_aggregate + curr_mu * curr_grad;
		curr_eta = (1 - curr_mu) * curr_epsilon + curr_mu * curr_alpha;
		curr_v = -1 * (sumsqr(curr_aggregate) + curr_eta);
		if curr_v > -precission
			break;
		end
		
		curr_y = curr_x - curr_aggregate;
		[dual_energy, curr_grad, upper_energy, labels_first] = dual(curr_y);

		if drawProfilePlot && any(iteration == profile_iterations)
			title = ['Profile aggregate bundle, iter = ', int2str(iteration)];
			if dual_energy <= curr_x_dual_energy + m * curr_v
				title = [title, ' (serious step)'];
			else
				title = [title, ' (null step)'];
			end
			profFig = profilePlot(@(point) -dual(point), curr_x, curr_y, 'title', title);
			profileCallback(profFig, iteration);
		end

		if dual_energy <= curr_x_dual_energy + m * curr_v
			curr_x = curr_y;
			curr_x_dual_energy = dual_energy;
			labels = labels_first;
		end
		curr_alpha = curr_x_dual_energy - dual_energy - curr_grad' * (curr_x - curr_y);
		curr_epsilon = curr_eta + curr_x_dual_energy - prev_x_dual_energy - curr_aggregate' * (curr_x - prev_x);

		prev_x = curr_x;
		prev_x_dual_energy = curr_x_dual_energy;
	end

	[~, energy, lowerBound, time] = wrapper.getState();
	oracle_calls = 1:(length(energy));
end

function [mu] = fundDirection(aggregate, grad, epsilon, alpha)
	p_norm_sqr = sumsqr(aggregate);
	g_norm_sqr = sumsqr(grad);
	p_times_g = aggregate' * grad;

	H = full((p_norm_sqr + g_norm_sqr - 2 * p_times_g));
	f = full((-p_norm_sqr + p_times_g - epsilon + alpha));
	lb = 0;
	ub = 1;
    % figure;
    % plot([0:0.1:1], arrayfun(@(mu) mu^2 * H / 2 + mu * f, [0:0.1:1]));
    precission = 1e-3;
    if abs(H) > precission
	    mu = - f / H;
	    if mu > ub
	    	mu = ub;
	    end
	    if mu < lb
	    	mu = lb;
	    end
    else
    	if f > 0
    		mu = lb;
    	else
    		mu = ub;
    	end
    end
end
