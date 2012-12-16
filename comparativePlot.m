function comparativePlot(name)
	% Compare different optimization methods and draw aggregate plot

	[unary, vertC, horC] = potentials(name);


	energy = cell(1);
	step = cell(1);
	lowerBound = cell(1);
	time = cell(1);
	legend_names = cell(1);
	plot_num = 1;

	% [labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
	% 																	struct('step', 0.01));
	% energy{plot_num} = curr_energy;
	% lowerBound{plot_num} = curr_lowerBound;
	% time{plot_num} = curr_time;
	% legend_names{plot_num} = 'Constant with step = 0.01';
	% plot_num = plot_num + 1;

	% [labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
	% 																	struct('step', 0.1));
	
	% energy{plot_num} = curr_energy;
	% lowerBound{plot_num} = curr_lowerBound;
	% time{plot_num} = curr_time;
	% legend_names{plot_num} = 'Constant with step = 0.1';
	% plot_num = plot_num + 1;

	% [labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
	% 																	struct('step', 1));
	% energy{plot_num} = curr_energy;
	% lowerBound{plot_num} = curr_lowerBound;
	% time{plot_num} = curr_time;
	% legend_names{plot_num} = 'Constant with step = 1';
	% plot_num = plot_num + 1;

	% [labels, curr_energy, curr_lowerBound, curr_time] = trwGridPotts(unary, vertC, horC, @constantSubgradient, ...
	% 																	struct('step', 3));
	% energy{plot_num} = curr_energy;
	% lowerBound{plot_num} = curr_lowerBound;
	% time{plot_num} = curr_time;
	% legend_names{plot_num} = 'Constant with step = 3';
	% plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time, curr_step] = trwGridPotts(unary, vertC, horC, ...
															@adaptiveSubgradient, struct());
	step{plot_num} = curr_step;
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	legend_names{plot_num} = 'Adaptive';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time, curr_step] = trwGridPotts(unary, vertC, horC, ...
															@backtrackingSubgradient, struct());
	step{plot_num} = curr_step;
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	legend_names{plot_num} = 'Backtracking';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time, curr_step] = trwGridPotts(unary, vertC, horC, ...
															@backtrackingSubgradient, struct('use_adaptive_init', 1));
	step{plot_num} = curr_step;
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	legend_names{plot_num} = 'Backtracking (from adaptive)';
	plot_num = plot_num + 1;

	[labels, curr_energy, curr_lowerBound, curr_time, curr_step] = trwGridPotts(unary, vertC, horC, ...
															@fletcherSubgradient, struct());
	step{plot_num} = curr_step;
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	legend_names{plot_num} = 'Fletcher';
	plot_num = plot_num + 1;

	% [labels, curr_energy, curr_time] = alphaExpansion(unary, vertC, horC);
	% energy{plot_num} = curr_energy;
	% lowerBound{plot_num} = curr_energy;
	% time{plot_num} = curr_time;
	% legend_names{plot_num} = 'Alpha expansion';
	% plot_num = plot_num + 1;


	[labels, curr_energy, curr_lowerBound, curr_time] = TRW_S(unary, vertC, horC);
	energy{plot_num} = curr_energy;
	lowerBound{plot_num} = curr_lowerBound;
	time{plot_num} = curr_time;
	legend_names{plot_num} = 'TRW-S';
	plot_num = plot_num + 1;
	% Set ylimits for small plot
	small_upper = curr_energy(5);
	small_lower = curr_lowerBound(5);




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


	fig_handle = figure;
	hold all;
	plot_arr = [];
	for i = 1:plot_num - 1
		curr_time = time{i};
		add = (curr_time(end) + 1) : max_time;
		add = reshape(add, length(add), 1);
		curr_time = [curr_time; add];
		curr_energy = [energy{i}; energy{i}(end) * ones(length(add), 1)];
		p = plot(curr_time, curr_energy);
		plotProperties(p);
		plot_arr = [plot_arr, p];
		curr_lowerBound = [lowerBound{i}; lowerBound{i}(end) * ones(length(add), 1)];
		p = plot(curr_time, curr_lowerBound, 'Color', get(p,'Color'));
		plotProperties(p);
	end
	

	hLegend = legend(plot_arr, legend_names);
	hTitle = title('Comparative plot');
	hXLabel = xlabel('Time (sec)');
	hYLabel = ylabel('Energy');
	set( gca                       , ...
	    'FontName'   , 'Helvetica' );
	set([hTitle, hXLabel, hYLabel], ...
	    'FontName'   , 'AvantGarde');
	set([hLegend]                 , ...
	    'FontSize'   , 11          );
	set([hXLabel, hYLabel]  , ...
	    'FontSize'   , 10          );
	set( hTitle                    , ...
	    'FontSize'   , 12          , ...
	    'FontWeight' , 'bold'      );
	hold off;

	out_filename = strcat('plots/comparative_', name);
	set(gcf, 'PaperPositionMode', 'auto');
	print('-depsc2', strcat(out_filename, '.eps'));
	saveas(fig_handle, out_filename,'fig');



	tmp = figure;
	hold all;
	plot_arr = [];
	for i = 1:plot_num - 1
		p = plot(energy{i});
		plotProperties(p);
		plot_arr = [plot_arr, p];
		p = plot(lowerBound{i}, 'Color', get(p,'Color'));
		plotProperties(p);
	end
	

	hLegend = legend(plot_arr, legend_names);
	hTitle = title('Comparative plot');
	hXLabel = xlabel('Outer iterations');
	hYLabel = ylabel('Energy');
	set( gca                       , ...
	    'FontName'   , 'Helvetica' );
	set([hTitle, hXLabel, hYLabel], ...
	    'FontName'   , 'AvantGarde');
	set([hLegend]                 , ...
	    'FontSize'   , 11          );
	set([hXLabel, hYLabel]  , ...
	    'FontSize'   , 10          );
	set( hTitle                    , ...
	    'FontSize'   , 12          , ...
	    'FontWeight' , 'bold'      );
	hold off;

	step_fig_handle = figure;
	hold all;
	plot_arr = [];
	for i = 1:plot_num - 1
		curr_time = time{i};
		add = (curr_time(end) + 1) : max_time;
		add = reshape(add, length(add), 1);
		curr_time = [curr_time; add];
		curr_step = [step{i}; step{i}(end) * ones(length(add), 1)];
		p = plot(curr_time, curr_step);
		plotProperties(p);
		plot_arr = [plot_arr, p];
	end
	

	hLegend = legend(plot_arr, legend_names);
	hTitle = title('Step plot');
	hXLabel = xlabel('Time (sec)');
	hYLabel = ylabel('Step size');
	set( gca                       , ...
	    'FontName'   , 'Helvetica' );
	set([hTitle, hXLabel, hYLabel], ...
	    'FontName'   , 'AvantGarde');
	set([hLegend]                 , ...
	    'FontSize'   , 11          );
	set([hXLabel, hYLabel]  , ...
	    'FontSize'   , 10          );
	set( hTitle                    , ...
	    'FontSize'   , 12          , ...
	    'FontWeight' , 'bold'      );
	hold off;


	out_filename = strcat('plots/step_', name);
	set(gcf, 'PaperPositionMode', 'auto');
	print('-depsc2', strcat(out_filename, '.eps'));
	saveas(step_fig_handle, out_filename,'fig');

	% out_filename = strcat(out_filename, '_small');
	% ylim([small_lower, small_upper]);
	% print('-depsc2', strcat(out_filename, '.eps'));
	% saveas(fig_handle, out_filename,'fig');
	% close;
end

function plotProperties(p)
	set(p, 'LineWidth', 1.5);
end
