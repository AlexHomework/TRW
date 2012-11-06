function [labels, energy, lowerBound, time, step] = dualDecomposition(K, N, f1, f2, dualStep, ...
													init_context, varargin)
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
	% Optional parameter lambda is for lambda_first initialization

	if length(varargin) > 0 && strcmp(varargin{1}, 'lambda') && ~isempty(varargin{2})
		lambda_first = varargin{2};
	else
		lambda_first = zeros(K, N);
	end
	lambda_second = -1 * lambda_first;
	lowerBound = [];
	energy = [];
	best_dual_energy = 0;
	time = [];
	step = [];
	context = init_context;
	t = cputime;
	for iteration = 1:10
		% Y minimization
		% The lower energy estimate
		[localEnergy, wholeEnergy, labels_first] = f1(lambda_first);
		dual_energy = localEnergy;
		upper_energy = wholeEnergy;
		[localEnergy, wholeEnergy, labels_second] = f2(lambda_second);
		dual_energy = dual_energy + localEnergy;
		upper_energy = min(wholeEnergy, upper_energy);

		best_dual_energy = max(best_dual_energy, dual_energy);
		lowerBound = [lowerBound, dual_energy];
		energy = [energy, upper_energy];


		[first_diff, context, curr_step] = dualStep(labels_first, ...
								labels_second, lowerBound, best_dual_energy, K, N, iteration, context);

		step = [step, curr_step];

		lambda_first = lambda_first + first_diff;
		lambda_second = lambda_second - first_diff;

		time = [time, cputime - t];
	end

	labels = labels_first;
end
