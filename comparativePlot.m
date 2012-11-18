function comparativePlot(unary, vertC, horC)
	% Compare different optimization methods and draw aggregate plot
	figure;
	hold all;
	plot_arr = [];
	legend_names = cell(1);
	plot_num = 1;

	[labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 0.01));
	X = [1:length(energy)];
	p = plot(X, energy);
	plot_arr = [plot_arr, p];
	plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 0.01';
	plot_num = plot_num + 1;

	[labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 0.1));
	X = [1:length(energy)];
	p = plot(X, energy);
	plot_arr = [plot_arr, p];
	plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 0.1';
	plot_num = plot_num + 1;

	[labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 1));
	X = [1:length(energy)];
	p = plot(X, energy);
	plot_arr = [plot_arr, p];
	plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 1';
	plot_num = plot_num + 1;

	[labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 3));
	X = [1:length(energy)];
	p = plot(X, energy);
	plot_arr = [plot_arr, p];
	plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 3';
	plot_num = plot_num + 1;

	[labels, energy, lowerBound, time] = trwGridPotts(unary, vertC, horC, @adaptiveSubgradient, struct());
	X = [1:length(energy)];
	p = plot(X, energy);
	plot_arr = [plot_arr, p];
	plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Adaptive';
	plot_num = plot_num + 1;

	[labels, energy, time] = alphaExpansion(unary, vertC, horC);
	plot_arr = [plot_arr, plot(X, double(energy) * ones(size(X)))];
	legend_names{plot_num} = 'Alpha expansion';
	plot_num = plot_num + 1;


	[labels, energy, lowerBound, time] = TRW_S(unary, vertC, horC);
	X = [1:length(energy)];
	p = plot(X, energy);
	plot_arr = [plot_arr, p];
	plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'TRW-S';
	plot_num = plot_num + 1;


	legend(plot_arr, legend_names);
	xlabel('Iteration')
	ylabel('Energy')
	hold off;
end