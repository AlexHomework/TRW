function compareEverything(experiment_name)
	% Compare different optimization methods, draw aggregate plots
	% and save everything (plots and the .mat file)
	% to experiments/experiment_name folder
	% 

	experiment_folder = strcat('experiments/', experiment_name);
	if exist(experiment_folder, 'dir') 
		disp('Experiment with such name already exist!');
		return;
	end

	mkdir(experiment_folder);

	data_set = getDataSet();
	for data_piece_ind = 1:length(data_set)
		data_piece = data_set{data_piece_ind};
		unary = data_piece.unary;
		vertC = data_piece.vertC;
		horC = data_piece.horC;
		data_name = data_piece.name;

		labels = cell(1);
		energy = cell(1);
		oracle_calls = cell(1);
		time = cell(1);
		step = cell(1);
		lowerBound = cell(1);
		legend_names = cell(1);
		plot_num = 1;

		baseline_algos = getBaselineAlgos();
		for algo_i = 1:length(baseline_algos)
			algo = baseline_algos{algo_i}{1};
			[curr_labels, curr_energy, curr_lowerBound, curr_time] = algo(unary, vertC, horC);
			labels{plot_num} = curr_labels;
			energy{plot_num} = curr_energy;
			lowerBound{plot_num} = curr_lowerBound;
			time{plot_num} = curr_time;
			step{plot_num} = [];
			oracle_calls{plot_num} = [];
			legend_names{plot_num} = baseline_algos{algo_i}{2};
			plot_num = plot_num + 1;
		end
		% Set ylimits for making small size plots
		small_upper = curr_energy(5);
		small_lower = curr_lowerBound(5);

		% Now test all TRW variants
		step_algos = getStepComputingAlgos();
		for algo_i = 1:length(step_algos)
			algo = step_algos{algo_i}{1};
			init_params = step_algos{algo_i}{2};
			[curr_labels, curr_energy, curr_lowerBound, ...
						curr_time, curr_step, curr_oracle_calls] = ...
						trwGridPotts(unary, vertC, horC, algo, init_params);
			labels{plot_num} = curr_labels;
			energy{plot_num} = curr_energy;
			lowerBound{plot_num} = curr_lowerBound;
			time{plot_num} = curr_time;
			step{plot_num} = curr_step;
			oracle_calls{plot_num} = curr_oracle_calls;
			legend_names{plot_num} = step_algos{algo_i}{3};
			plot_num = plot_num + 1;
		end

		plot_num = plot_num - 1;





		max_time = 0;
		max_iter = 0;
		for i = 1:plot_num
			curr_time = time{i};
			if curr_time(end) > max_time
				max_time = curr_time(end);
			end
			if length(curr_time) > max_iter
				max_iter = length(curr_time);
			end
		end

		time_fig = figure;
		hold all;
		plot_arr = [];
		curr_legend_names = {};
		for i = 1:plot_num
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
			curr_legend_names{end + 1} = legend_names{i};
			plotProperties(p);
		end
		figureProperties(plot_arr, curr_legend_names, 'Time plot', 'Time (sec)', 'Energy');
		hold off;

		out_filename = strcat(experiment_folder, '/time_', data_name);
		set(gcf, 'PaperPositionMode', 'auto');
		print('-depsc2', strcat(out_filename, '.eps'));
		saveas(time_fig, out_filename,'fig');

		out_filename = strcat(out_filename, '_small');
		ylim([small_lower, small_upper]);
		print('-depsc2', strcat(out_filename, '.eps'));
		saveas(time_fig, out_filename,'fig');
		close;



		outer_iterations_fig = figure;
		hold all;
		plot_arr = [];
		curr_legend_names = {};
		for i = 1:plot_num
			if (length(energy{i}) > 1)
				p = plot(energy{i});
				plotProperties(p);
				plot_arr = [plot_arr, p];
				p = plot(lowerBound{i}, 'Color', get(p,'Color'));
				curr_legend_names{end + 1} = legend_names{i};
				plotProperties(p);
			end
		end
		figureProperties(plot_arr, curr_legend_names, 'Comparative plot', 'Outer iterations', 'Energy');
		hold off;

		out_filename = strcat(experiment_folder, '/outer_iterations_', data_name);
		set(gcf, 'PaperPositionMode', 'auto');
		print('-depsc2', strcat(out_filename, '.eps'));
		saveas(outer_iterations_fig, out_filename,'fig');

		out_filename = strcat(out_filename, '_small');
		ylim([small_lower, small_upper]);
		print('-depsc2', strcat(out_filename, '.eps'));
		saveas(outer_iterations_fig, out_filename,'fig');
		close;


		step_fig = figure;
		hold all;
		plot_arr = [];
		curr_legend_names = {};
		for i = 1:plot_num
			if (length(step{i}) > 0)
				curr_time = time{i};
				add = (curr_time(end) + 1) : max_time;
				add = reshape(add, length(add), 1);
				curr_time = [curr_time; add];
				curr_step = [step{i}; step{i}(end) * ones(length(add), 1)];
				p = plot(curr_time, curr_step);
				curr_legend_names{end + 1} = legend_names{i};
				plotProperties(p);
				plot_arr = [plot_arr, p];
			end
		end
		figureProperties(plot_arr, curr_legend_names, 'Step plot', 'Step size', 'Energy');
		hold off;


		out_filename = strcat(experiment_folder, '/step_', data_name);
		set(gcf, 'PaperPositionMode', 'auto');
		print('-depsc2', strcat(out_filename, '.eps'));
		saveas(step_fig, out_filename,'fig');
		close;


		

		oracle_fig = figure;
		hold all;
		plot_arr = [];
		curr_legend_names = {};
		for i = 1:plot_num
			if (length(oracle_calls{i}) > 0)
				curr_oracle_calls = oracle_calls{i};
				curr_energy = energy{i};
				p = plot(curr_oracle_calls, curr_energy);
				plotProperties(p);
				plot_arr = [plot_arr, p];
				curr_lowerBound = lowerBound{i};
				p = plot(curr_oracle_calls, curr_lowerBound, 'Color', get(p,'Color'));
				curr_legend_names{end + 1} = legend_names{i};
				plotProperties(p);
			end
		end
		figureProperties(plot_arr, curr_legend_names, 'Oracle calls plot', 'Oracle calls', 'Energy');
		hold off;

		out_filename = strcat(experiment_folder, '/oracle_', data_name);
		set(gcf, 'PaperPositionMode', 'auto');
		print('-depsc2', strcat(out_filename, '.eps'));
		saveas(oracle_fig, out_filename,'fig');

		out_filename = strcat(out_filename, '_small');
		ylim([small_lower, small_upper]);
		print('-depsc2', strcat(out_filename, '.eps'));
		saveas(oracle_fig, out_filename,'fig');
		close;


		out_filename = strcat(experiment_folder, '/data.mat');
		save(out_filename, 'labels', 'energy', 'oracle_calls', 'time', 'step', 'lowerBound', 'legend_names');
	end
end

function plotProperties(p)
	set(p, 'LineWidth', 1.5);
end

function figureProperties(plot_arr, legend_names, title_str, xlabel_str, ylabel_str)
	hLegend = legend(plot_arr, legend_names);
	hTitle = title(title_str);
	hXLabel = xlabel(xlabel_str);
	hYLabel = ylabel(ylabel_str);
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
end
