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

	function [labels, energy, lowerBound, time] = alphaExpansionWithLower(unary, vertC, horC)
		% Add dummy lower bound to alpha expansion
		[labels, energy, time] = alphaExpansion(unary, vertC, horC);
		lowerBound = energy;
	end

	algorithms = {
		{@alphaExpansionWithLower, 'Alpha Expansion'},...
		{@TRW_S, 'TRW-S'},...
	};

end
