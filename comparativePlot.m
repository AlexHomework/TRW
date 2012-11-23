function comparativePlot(unary, vertC, horC)
	% Compare different optimization methods and draw aggregate plot
	energy = cell(1);
	lowerBound = cell(1);
	time = cell(1);
	legend_names = cell(1);
	plot_num = 1;

	[labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 0.01));
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	% X = [1:length(energy)];
	% p = plot(X, energy);
	% plot_arr = [plot_arr, p];
	% plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 0.01';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 0.1));
	
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	% X = [1:length(energy)];
	% p = plot(X, energy);
	% plot_arr = [plot_arr, p];
	% plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 0.1';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 1));
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	% X = [1:length(energy)];
	% p = plot(X, energy);
	% plot_arr = [plot_arr, p];
	% plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 1';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
																		struct('step', 3));
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	% X = [1:length(energy)];
	% p = plot(X, energy);
	% plot_arr = [plot_arr, p];
	% plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Constant with step = 3';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @adaptiveSubgradient, struct());
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	% X = [1:length(energy)];
	% p = plot(X, energy);
	% plot_arr = [plot_arr, p];
	% plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'Adaptive';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_time] = alphaExpansion(unary, vertC, horC);
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_energy;
	time{plot_num} = curr_time;
	% plot_arr = [plot_arr, plot(X, double(energy) * ones(size(X)))];
	legend_names{plot_num} = 'Alpha expansion';
	plot_num = plot_num + 1;


	[labels, curr_energy, curr_lowerBound, curr_time] = TRW_S(unary, vertC, horC);
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	% X = [1:length(energy)];
	% p = plot(X, energy);
	% plot_arr = [plot_arr, p];
	% plot(X, lowerBound, 'Color', get(p,'Color'));
	legend_names{plot_num} = 'TRW-S';
	plot_num = plot_num + 1;




	max_time = 0;
	max_iter = 0;
	for i = 1:plot_num - 1
		curr_time = time{i};
		if curr_time(end) > max_time
			max_time = curr_time(end);
		end
		if length(curr_time) > max_iter
			max_iter = length(curr_time);
		end
	end


	figure;
	hold all;
	plot_arr = [];
	for i = 1:plot_num - 1
		curr_time = time{i};
		add = (curr_time(end) + 1) : max_time;
		add = reshape(add, length(add), 1);
		curr_time = [curr_time; add];
		curr_energy = [energy{i}; energy{i}(end) * ones(length(add), 1)];
		p = plot(curr_time, curr_energy);
		plot_arr = [plot_arr, p];
		curr_lowerBound = [lowerBound{i}; lowerBound{i}(end) * ones(length(add), 1)];
		plot(curr_time, curr_lowerBound, 'Color', get(p,'Color'));
	end
	
	legend(plot_arr, legend_names);
	xlabel('Time')
	ylabel('Energy')
	hold off;

	% % Iteration/energy plot
	% figure;
	% hold all;
	% plot_arr = [];
	% for i = 1:plot_num - 1
	% 	curr_X = [1:length(energy{i})];
	% 	add = (length(energy{i}) + 1) : max_iter;
	% 	curr_X = [curr_X, add];
	% 	curr_X = reshape(curr_X, length(curr_X), 1);
	% 	curr_energy = [energy{i}; energy{i}(end) * ones(length(add), 1)];
	% 	p = plot(curr_X, curr_energy);
	% 	plot_arr = [plot_arr, p];
	% 	curr_lowerBound = [lowerBound{i}; lowerBound{i}(end) * ones(length(add), 1)];
	% 	plot(curr_X, curr_lowerBound, 'Color', get(p,'Color'));
	% end

	% legend(plot_arr, legend_names);
	% xlabel('Iteration')
	% ylabel('Energy')
	% hold off;
end