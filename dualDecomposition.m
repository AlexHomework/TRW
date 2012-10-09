function [labels, energy, lowerBound, time] = dualDecomposition(K, N, f1, f2)
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

	gamma0 = 1.5;
	gamma1 = 0.5;
	epsilon = @(n) 1 / n;
	delta_prev = 1000;

	lambda_first = zeros(K, N);
	lambda_second = zeros(K, N);
	lowerBound = [];
	energy = [];
	best_dual_energy = 0;
	time = [];
	t = cputime;
	for iteration = 1:1
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


		 % Adaptive projected subgradient step computation
		if iteration == 1
			delta = delta_prev;
		else
			if dual_energy > lowerBound(iteration - 1)
				delta = gamma0 * delta_prev;
			else
				delta = max(gamma1 * delta_prev, epsilon(iteration));
			end
		end
		delta_prev = delta;
		alpha_n = best_dual_energy + delta - dual_energy;
		alpha_n = alpha_n / sum(sum(labels_first ~= labels_second));


		% Lambda projected subgradient maximization
		for p = 1:K
			lambda_first(p, :) = lambda_first(p, :) + alpha_n * ((labels_first == p) - (labels_second == p));
			lambda_second(p, :) = lambda_second(p, :) + alpha_n * ((labels_second == p) - (labels_first == p));
		end

		time = [time, cputime - t];
	end

	labels = labels_first;
end
