function [labels, energy, lowerBound, time] = bundleDual(unary, vertC, horC, varargin)
	% Optimize dual energy via bundle method
	% 

	[K, N, M] = size(unary);
	[epsilon, mL, gamma, wMax, maxBundleSize, iterations_count] = process_options(varargin, ...
											'eps', 0.001, 'mL', 0.1, 'gamma', 0.1, 'wMax', 10, ...
											'bundleSize', 10, 'iter', 700);

	lambda_first = zeros(K * N * M, 1);
	
	wrapper = gridDualWrapper(unary, vertC, horC);
	function [dual_energy, grad, upper_energy, labels_first, labels_second] = dual(lambda)
		[dual_energy, grad, upper_energy, labels_first, labels_second] = ...
								wrapper.dual(reshape(lambda, K, N * M));
		grad = grad(:);
	end


	clear bundle;
	[dual_energy, grad, upper_energy, labels_first, labels_second] = dual(lambda_first);
	min_upper_energy = upper_energy;
	max_dual_energy = dual_energy;
	bundle.f = [dual_energy];
	bundle.g = grad;
	bundle.dotProduct = grad' * lambda_first;
	bundle.gNorm = sumsqr(grad);
	bundle.size = 1;
	lambdaCenter = lambda_first;
	f_center = dual_energy;
	w = 1;
	for iteration = 1:iterations_count
		% Firstly get current bundle maximization
		% point and value of bundle funcion in this point.
		[lambda_next, bundle_val_lambda_next, bundle] = maximizeBundle(bundle, lambdaCenter, w, maxBundleSize);
		delta = bundle_val_lambda_next - f_center;
		[dual_energy, grad, upper_energy, labels_first, labels_second] = dual(lambda_next);
		min_upper_energy = min([min_upper_energy, upper_energy]);
		max_dual_energy = max([max_dual_energy, dual_energy]);
		if ((dual_energy - f_center) >= mL * (bundle_val_lambda_next - f_center))
			lambdaCenter = lambda_next;
			f_center = dual_energy;
			w = gamma * (min_upper_energy - max_dual_energy) /  norm(grad);
			if w > wMax
				w = wMax;
			end
		end


		% Update bundle
		bundle.f(end + 1, 1) = dual_energy;
		bundle.g(:, end + 1) = grad;
		bundle.dotProduct(end + 1, 1) = grad' * lambda_next;
		bundle.gNorm(end + 1, 1) = sumsqr(grad);
		bundle.size = bundle.size + 1;


		if delta < epsilon
			sprintf('Delta is less then epsilon (''%d'' < ''%d''), finishing on ''%d'' iteration.', delta, epsilon, iteration)
			break
		end

	end

	labels = labels_first;
	[~, energy, lowerBound, time] = wrapper.getState();
end

function [lambdaMax, value, cleanedBundle] = maximizeBundle(bundle, lambdaCenter, w, maxBundleSize)
	% Solve bundle approximation maximization problem
	% and possible remove least violated constraints (if
	% bundle exceeded maximum size).
	% 
	% bundle.f(i) is the function value on the i-th step (in the lambda(:, i) point)
	% bundle.g(:, i) is the gradient vector on the i-th step
	% bundle.gNorm(i) = sumsqr(bundle.g(:, i))
	% bundle.dotProduct(i) is bundle.g(:, i)' * [point on the i-th step]
	% lambdaCenter is the regularization center
	% w is regularization constant
	% 

	N = bundle.size;
	nFeatures = length(lambdaCenter);
	H = diag(bundle.gNorm) / w;
	hyperplane = -bundle.dotProduct' +  lambdaCenter' * bundle.g;
	f = hyperplane(:) + bundle.f;

	Aeq = ones(1, N);
	beq = 1;
	lb = zeros(N, 1);
	ub = ones(N, 1);

    options = optimset('Algorithm', 'interior-point-convex');
    xi = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);


    lambdaDiff = bundle.g * xi(:) / w;
    lambdaMax = lambdaDiff + lambdaCenter;
    constraints = hyperplane + lambdaDiff' * bundle.g;
    constraints = constraints(:) + bundle.f;
    if bundle.size >= maxBundleSize
    	% Remove least violated constraint
    	clear cleanedBundle;
    	[~, worst] = max(constraints);
    	idx = [1:(worst - 1), (worst + 1):N];
    	cleanedBundle.f = bundle.f(idx);
    	cleanedBundle.g = bundle.g(:, idx);
    	cleanedBundle.dotProduct = bundle.dotProduct(idx);
    	cleanedBundle.gNorm = bundle.gNorm(idx);
    	cleanedBundle.size = bundle.size - 1;
    else
	    cleanedBundle = bundle;
	end
    value = min(constraints);
end


