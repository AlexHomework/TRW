function [labels, primal_energy, dual_energy, time, oracle_calls] = bfgsDual(unary, vertC, horC, varargin)
	addpath('hanso/');

	[max_iter, drawProfilePlot, profile_iterations, profileCallback] = process_options(varargin, 'maxIter', 30, 'profile', false, ...
																	  'profileIters', [5, 10], 'profileCallback', @(fig, iter) 0);

	if drawProfilePlot
		disp(['Profile plots will be generated. Be careful, it', ...
			  ' can take a while and thus distort time measurements.']);
	end
	[K, N, M] = size(unary);
	pars.nvar = K * N * M;
	wrapper = gridDualWrapper(unary, vertC, horC, 'save_points', profile_iterations);
	profile_from = [];
	profile_to = [];
	function [value, derivative] = minus_dual(lambda, pars)
		[value, derivative] = wrapper.dual(reshape(lambda, K, N * M));
		value = -value;
		derivative = -derivative(:);
	end
	pars.fgname = @minus_dual;
	options.prtlevel = 2;
	options.maxit = max_iter;
	options.x0 = zeros(K * N * M, 1);
	lambda = hanso(pars, options);
	if ~drawProfilePlot
		[~, primal_energy, dual_energy, time] = wrapper.getState();
	else
		[~, primal_energy, dual_energy, time, lambdaHist] = wrapper.getState();

		for i = 1:length(profile_iterations)
			curr = profile_iterations(i);
			if (curr > numel(lambdaHist))
				break
			end
			title = ['Profile bfgs, iter = ', int2str(curr)];
			profFig = profilePlot(@wrapper.dual, lambdaHist{curr}, lambdaHist{curr + 1}, 'title', title);
			profileCallback(profFig, curr);
		end
		clear lambdaHist;
	end
	[~, ~, ~, labels] = gridDual(reshape(lambda, K, N * M), unary, vertC, horC);
	oracle_calls = 1:length(dual_energy);
end
