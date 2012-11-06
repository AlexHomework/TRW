function stepPlot(unary, vertC, horC, dualStep, init_context)
	% Plot dual step size versus iteration
	figure;
	[labels, energy, lowerBound, time, step] = trwGridPotts(unary, vertC, horC, dualStep, init_context);
	X = [1:length(energy)];
	plot(X, step);
	xlabel('Iteration');
	ylabel('Step size');
end