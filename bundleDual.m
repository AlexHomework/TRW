function [labels, energy, lowerBound, time, oracle_calls] = bundleDual(unary, vertC, horC, varargin)
	% Optimize dual energy via bundle method
	% 

	[K, N, M] = size(unary);
	[epsilon, mL, gamma, wMax, maxBundleSize, iterations_count, drawProfilePlot, profile_iterations, profileCallback] = ...
											process_options(varargin, ...
											'eps', 0.001, 'mL', 0.1, 'gamma', 0.1, 'wMax', 10, ...
											'bundleSize', 10, 'maxIter', 40, 'profile', false, ...
											'profileIters', [5, 10], 'profileCallback', @(fig, iter) 0);


    wMin = 1e-10;
	lambda_first = zeros(K * N * M, 1);
	
	if drawProfilePlot
		disp(['Profile plots will be generated. Be careful, it', ...
			  ' can take a while and thus distort time measurements.']);
	end

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
	bundle.size = 1;
	lambdaCenter = lambda_first;
	f_center = dual_energy;
	w = 1;
	for iteration = 1:iterations_count
		% Firstly get current bundle maximization
		% point and value of bundle funcion at this point.
		[lambda_next, bundle_val_lambda_next, bundle] = maximizeBundle(bundle, lambdaCenter, w, maxBundleSize);
		
		if drawProfilePlot && any(iteration == profile_iterations)
			direction = lambda_next - lambdaCenter;
			title = ['Profile (bundle size = ', int2str(maxBundleSize)];
			title = [title, ', iter = ', int2str(iteration), ', w = ', num2str(w), ')'];
			profFig = profilePlot(@dual, lambdaCenter, lambda_next, 'title', title);
			moments = [0:0.1:2];
			for bundleIdx = 1:bundle.size
				bundleValues = arrayfun(@(x) bundle.f(bundleIdx) + bundle.g(:, bundleIdx)' *(lambdaCenter + x * direction) - bundle.dotProduct(bundleIdx), moments);
				p3 = plot(moments, bundleValues);
			end
			f = bundle.f - bundle.dotProduct + bundle.g' * lambdaCenter;
			direction_proj = bundle.g' * direction;
			values = arrayfun(@(x) min(f + x * direction_proj), moments);
			p2 = plot(moments, values, '-g');
			profileCallback(profFig, iteration);
		end

		delta = bundle_val_lambda_next - f_center;
		[dual_energy, grad, upper_energy, labels_first, labels_second] = dual(lambda_next);
		min_upper_energy = min([min_upper_energy, upper_energy]);
		max_dual_energy = max([max_dual_energy, dual_energy]);
		if ((dual_energy - f_center) >= mL * (bundle_val_lambda_next - f_center))
			lambdaCenter = lambda_next;
			f_center = dual_energy;
			w = norm(grad) / (gamma * (min_upper_energy - max_dual_energy));
			if w > wMax
				w = wMax;
			end
			if w < wMin
				w = wMin;
			end
		end


		% Update bundle
		bundle.f(end + 1, 1) = dual_energy;
		bundle.g(:, end + 1) = grad;
		bundle.dotProduct(end + 1, 1) = grad' * lambda_next;
		bundle.size = bundle.size + 1;

		


		if delta < epsilon
			sprintf('Delta is less then epsilon (''%d'' < ''%d''), finishing on ''%d'' iteration.', delta, epsilon, iteration)
			break
		end

	end

	labels = labels_first;
	[~, energy, lowerBound, time] = wrapper.getState();
	oracle_calls = 1:(length(energy));
end

function [lambdaMax, value, cleanedBundle] = maximizeBundle(bundle, lambdaCenter, w, maxBundleSize)
	% Solve bundle approximation maximization problem
	% and possible remove least violated constraints (if
	% bundle exceeded maximum size).
	% 
	% bundle.f(i) is the function value on the i-th step (in the lambda(:, i) point)
	% bundle.g(:, i) is the gradient vector on the i-th step
	% bundle.dotProduct(i) is bundle.g(:, i)' * [point on the i-th step]
	% lambdaCenter is the regularization center
	% w is regularization constant
	% 

	N = bundle.size;
	nFeatures = length(lambdaCenter);
	H = bundle.g' * bundle.g / w;
	f = bundle.f - bundle.dotProduct + bundle.g' * lambdaCenter;

	Aeq = ones(1, N);
	beq = 1;
	lb = zeros(N, 1);
	ub = ones(N, 1);

    options = optimset('Algorithm', 'interior-point-convex', 'Display', 'off');
    xi = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);


    lambdaDiff = bundle.g * xi(:) / w;
    lambdaMax = lambdaDiff + lambdaCenter;
    constraints = f + bundle.g' * lambdaDiff;
    value = min(constraints);

    % % Compare fast optimization with the slow one to make sure that everything is implemented correctly.
    % [lambdaMaxOld, valueOld] = maximizeBundleOld(bundle, lambdaCenter, w, maxBundleSize);
    % diff = norm(lambdaMaxOld-lambdaMax)/norm(lambdaMaxOld+lambdaMax);
    % disp(diff);
    % disp(value)
    % disp(valueOld)

    if bundle.size >= maxBundleSize
    	% Remove least violated constraint
    	clear cleanedBundle;
    	[~, worst] = max(constraints);
    	idx = [1:(worst - 1), (worst + 1):N];
    	cleanedBundle.f = bundle.f(idx);
    	cleanedBundle.g = bundle.g(:, idx);
    	cleanedBundle.dotProduct = bundle.dotProduct(idx);
    	cleanedBundle.size = bundle.size - 1;
    else
	    cleanedBundle = bundle;
	end
end


function [lambdaMax, value] = maximizeBundleOld(bundle, lambdaCenter, w, maxBundleSize)
	N = bundle.size;
	nFeatures = length(lambdaCenter);
	H = w * speye(nFeatures + 1);
	H(nFeatures + 1, nFeatures + 1) = 0;
	f = -w * lambdaCenter;
	f(end + 1) = -1;
	A = [-bundle.g', ones(N, 1)];
	b = -bundle.dotProduct + bundle.f;
	options = optimset('Algorithm', 'interior-point-convex');
	vars = quadprog(H, f, A, b, [], [], [], [], [], options);
	lambdaMax = vars(1:end - 1);
	value = vars(end);
end
