function [labels, energy, lowerBound, time] = bundleDual(unary, vertC, horC, varargin)

	[K, N, M] = size(unary);


	epsilon = 0.001;
	mL = 0.1;
	gamma = 0.1;
	wMax = 10;
	maxBundleSize = 10;

	iterations_count = 200;

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
	bundle.lambda = lambda_first;
	bundle.f = [dual_energy];
	bundle.g = grad;
	bundle.size = 1;
	lambda_center = lambda_first;
	f_center = dual_energy;
	w = 1;
	for iteration = 1:iterations_count
		% Firstly get current bundle maximization
		% point and value of bundle funcion in this point.
		[lambda_next, bundle_val_lambda_next, bundle] = maximizeBundle(bundle, lambda_center, w, maxBundleSize);
		delta = bundle_val_lambda_next - f_center;
		[dual_energy, grad, upper_energy, labels_first, labels_second] = dual(lambda_next);
		min_upper_energy = min([min_upper_energy, upper_energy]);
		max_dual_energy = max([max_dual_energy, dual_energy]);
		if ((dual_energy - f_center) >= mL * (bundle_val_lambda_next - f_center))
			lambda_center = lambda_next;
			f_center = dual_energy;
			w = gamma * (min_upper_energy - max_dual_energy) /  norm(grad);
			if w > wMax
				w = wMax;
			end
		end


		% Update bundle
		bundle.lambda(:, end + 1) = lambda_next;
		bundle.f(end + 1, 1) = dual_energy;
		bundle.g(:, end + 1) = grad;
		bundle.size = bundle.size + 1;


		if delta < epsilon
			disp('delta is less then epsilon');
			delta
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
	% bundle.lambda(:, i) is point on the i-th step
	% bundle.f(i) is the function value on the i-th step (in the lambda(:, i) point)
	% bundle.g(:, i) is the gradient vector on the i-th step
	% lambdaCenter is the regularization center
	% w is regularization constant
	% 

	N = bundle.size;

	% nFeatures = length(lambdaCenter);
	% H = (w / 2) * speye(nFeatures + 1);
	% H(nFeatures + 1, nFeatures + 1) = 0;
	% f = -w * lambdaCenter;
	% f(end + 1) = -1;
	% A = [-bundle.g', ones(N, 1)];
	% hyperplane = -sum(bundle.g .* bundle.lambda);
	% b = hyperplane(:) + bundle.f;
 %    options = optimset('Algorithm', 'interior-point-convex');
 %    vars = quadprog(H, f, A, b, [], [], [], [], [], options);
 %    lambda = vars(1:end - 1)';
 %    lambdaMax = vars(end);



	gNorms = sum(bundle.g.^2);
	H = diag(gNorms) / w;
	hyperplane = sum(bundle.g .* (bsxfun(@plus, -bundle.lambda, lambdaCenter)));
	f = hyperplane(:) + bundle.f;

	Aeq = ones(1, N);
	beq = 1;
	lb = zeros(N, 1);
	ub = ones(N, 1);

    options = optimset('Algorithm', 'interior-point-convex');
    xi = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);


    lambdaMax = bundle.g * xi(:) / w + lambdaCenter;
    constraints = sum(bundle.g .* (bsxfun(@plus, -bundle.lambda, lambdaMax)));
    constraints = constraints(:) + bundle.f;
    while bundle.size >= maxBundleSize
    	% Remove least violated constraint
    	[~, worst] = max(constraints);
        constraints(worst) = [];
    	bundle.lambda(:, worst) = [];
    	bundle.f(worst) = [];
    	bundle.g(:, worst) = [];
    	bundle.size = bundle.size - 1;
    end
    cleanedBundle = bundle;
    value = min(constraints);
end


