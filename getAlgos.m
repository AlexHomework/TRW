function algorithms = getBaselineAlgos()
	% Return cell array of baseline optimization functions
	% which we will compare to ours.
	% 
	% algorithms{i}{1} is the i-th function
	% algorithms{i}{2} is the i-th function name (for plots)
	% 
	% Each function should take
	% (unary, vertC, horC) arguments and return
	% [labels, energy, lowerBound, time]
	%

	function [labels, energy, lowerBound, time, oracle_calls] = alphaExpansionWithLower(unary, vertC, horC)
		% Add dummy lower bound and oracle calls arrays to alpha expansion
		[labels, energy, time] = alphaExpansion(unary, vertC, horC);
		lowerBound = energy;
		oracle_calls = [];
	end

	function [labels, energy, lowerBound, time, oracle_calls] = TRW_SWithOracle(unary, vertC, horC)
		% Add dummy oracle calls array to TRW_S
		[labels, energy, lowerBound, time] = TRW_S(unary, vertC, horC);
		oracle_calls = [];
	end

	algorithms = {
		{@bfgsDual, 'bfgs'},...
		{@(u, v, h) bundleDual(u, v, h, 'bundleSize', 10), 'Bundle (10)'},...
		{@(u, v, h) bundleDual(u, v, h, 'bundleSize', 10), 'Bundle (100)'},...
		{@alphaExpansionWithLower, 'Alpha Expansion'},...
		{@TRW_SWithOracle, 'TRW-S'},...
	};

end
