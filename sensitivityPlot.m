function sensitivityPlot(unary, vertC, horC, dualStep, init_context)
	% Plot energy versus iteration for different random initial lambda

	% Count of random start runs
	N = 5;

	figure;
	hold on;

	[labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, dualStep, init_context);
	X = [1:length(energy)];
	p1 = plot(X, energy, '-b');
	plot(X, lowerBound, '-b');
	for i = 1:N
		[labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, dualStep, ...
																		init_context, 'random_init');
		p2 = plot(X, energy, '-r');
		plot(X, lowerBound, '-r');
	end
	legend([p1, p2], 'Zero start', 'Random start');

	xlabel('Iteration');
	ylabel('Energy');
	hold off;
end