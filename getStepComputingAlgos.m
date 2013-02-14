function algorithms = getStepComputingAlgos()
	% Return cell array of functions to try in
	% 1-d step optimization of TRW dual function optimization.
	% 
	% algorithms{i}{1} is the i-th function
	% algorithms{i}{2} is the i-th function init parameters
	% algorithms{i}{3} is the i-th function name (for plots)
	% 
	% Each function should take
	% (func, grad, lower_bound, iteration, context) arguments and return
	% [context, alpha_n, f_n]
	% Where func  -- is our 1-d function (TRW dual in optimization direction)
	% grad        -- subgradient at the base point ( func(0) )
	% lower_bound -- values of our target function (TRW dual)
	%		on the previous iterations
	% iteration   -- current iteration number
	% context     -- space where function can store intermediate computations between iterations
	%
	% And alpha_n is chosen step, f_n = func(alpha_n)
	% 

	function [labels, energy, lowerBound, time] = alphaExpansionWithLower(unary, vertC, horC)
		% Add dummy lower bound to alpha expansion
		[labels, energy, time] = alphaExpansion(unary, vertC, horC);
		lowerBound = energy;
	end

	algorithms = {
		% {@constantSubgradient, struct('step', 0.01), 'Constant 0.01'}, ...
		% {@constantSubgradient, struct('step', 0.1), 'Constant 0.1'}, ...
		% {@constantSubgradient, struct('step', 0.1), 'Constant 0.1'}, ...
		% {@constantSubgradient, struct('step', 1), 'Constant 1'}, ...
		% {@constantSubgradient, struct('step', 3), 'Constant 3'}, ...
		{@adaptiveSubgradient, struct(), 'Adaptive'}, ...
		{@backtrackingSubgradient, struct(), 'Backtracking'}, ...
		{@backtrackingSubgradient, struct('use_adaptive_init', 1), 'Backtracking (from adative)'}, ...
		{@fletcherSubgradient, struct(), 'Fletcher'}, ...
	};

end
