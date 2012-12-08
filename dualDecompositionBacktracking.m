function [labels, energy, lowerBound, time, step] = dualDecompositionBacktracking(K, N, f1, f2, varargin)
	% Equivalent to dualDecomposition.m but
	% using Backtracking method to optimize dual step

	% Whole dual problem is to find
	% min f1'(x1, \theta1) + min f2'(x2, \theta2)
	% subject to x1 = x2
	% f1(\lambda) = min over x1 (f1'(x1, \theta1 + \lambda))
	% 
	% \theta is a matrix:
	% [K, N] = size(\theta)
	% 
	% [localEnergy, wholeEnergy, labels] = f1(\lambda)
	% where labels = argmin f1'(x, \theta1 + \lambda)
	% localEnergy  = min f1'(x, \theta1 + \lambda)
	% wholeEnergy  = f1'(labels, \theta1 + \lambda) + f2'(labels, \theta1 + \lambda)
	%
	% Optional parameter lambda is for lambda initialization

	function [dual_energy, upper_energy, labels_first, labels_second] = minimizeDual(lambda)
		[localEnergy, wholeEnergy, labels_first] = f1(lambda);
		dual_energy = localEnergy;
		upper_energy = wholeEnergy;
		[localEnergy, wholeEnergy, labels_second] = f2(-1 * lambda);
		dual_energy = dual_energy + localEnergy;
		upper_energy = min(wholeEnergy, upper_energy);
	end
		

	if ~isempty(varargin) & strcmp(varargin{1}, 'lambda') & ~isempty(varargin{2})
		lambda = varargin{2};
	else
		lambda = zeros(K, N);
	end
	t = cputime;


	% The lower energy estimate
	[dual_energy, upper_energy, labels_first, labels_second] = minimizeDual(lambda);
	best_dual_energy = dual_energy;
	lowerBound = [dual_energy];
	energy = [upper_energy];
	time = [cputime - t];
	step = [0];
	for iteration = 1:100
		grad = zeros(K, N);
		for p = 1:K
			grad(p, :) = reshape(((labels_first == p) - (labels_second == p)), 1, N);
		end

		[curr_step, dual_energy, upper_energy, labels_first, labels_second] = ...
						minBacktracking(@(step) minimizeDual(lambda + step * grad), dual_energy, ...
																		sum(sum(abs(grad))), 1);

		best_dual_energy = max(best_dual_energy, dual_energy);
		lowerBound = [lowerBound; dual_energy];
		energy = [energy; upper_energy];

		step = [step; curr_step];

		lambda = lambda + curr_step * grad;

		time = [time; cputime - t];
	end

	labels = labels_first;
end

function [best_h, dual_energy, upper_energy, labels_first, labels_second] = ...
													minBacktracking(func, func_0, grad_sq, h_init)
	max_iter = 15;
	fact = 0.7;
	rho = 0.1;
	h = h_init;
	it = 1;
	[dual_energy, upper_energy, labels_first, labels_second] = func(h);

	% dual_energy
	% func_0
	% grad_sq
	% h
	% dual_energy - (func_0 + rho * h * grad_sq)
	% pause
	while (it <= max_iter) & (dual_energy < (func_0 + rho * h * grad_sq))
		[dual_energy, upper_energy, labels_first, labels_second] = func(h);
		h = h * fact;
		it = it + 1;
	end
	best_h = h;
	if it > max_iter
		disp('Reached maximum number of iterations!');
	else
		disp('Approximate optimization succeeded.');
	end
	% dual_energy
	% func_0
	% grad_sq
	% h
end