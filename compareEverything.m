function compareEverything(experiment_name)
	% Compare different optimization methods, draw aggregate plots
	% and save everything (plots and the .mat file)
	% to experiments/experiment_name folder
	% 

	experiment_folder = strcat('experiments/', experiment_name);

	data_set = getDataSet();
	baseline_algos = getBaselineAlgos();
	step_algos = getStepComputingAlgos();
	if exist(experiment_folder, 'dir')
		% Find that we have already computed
		processed_filename = strcat(experiment_folder, '/processed_list.mat');
		load(processed_filename, 'data_set_names', 'baseline_algos_names', 'step_algos_names');
	else
		data_set_names = {};
		baseline_algos_names = {};
		step_algos_names = {};
		mkdir(experiment_folder);
	end

	for data_piece_ind = 1:length(data_set)
		data_piece = data_set{data_piece_ind};
		unary = data_piece.unary;
		vertC = data_piece.vertC;
		horC = data_piece.horC;
		data_name = data_piece.name;
		% known_data_piece == true iff current data piece
		% was already processed with some algorithms.
		known_data_piece = sum(ismember(data_set_names, data_name)) == 1;

		if (known_data_piece)
			results_filename = strcat(experiment_folder, '/data_', data_name, '.mat');
			load(results_filename, 'labels', 'energy', 'oracle_calls', 'time', 'step', 'lowerBound', 'legend_names');
			plot_num = length(labels) + 1;
		else
			% Start everything from scratch.
			labels = cell(1);
			energy = cell(1);
			oracle_calls = cell(1);
			time = cell(1);
			step = cell(1);
			lowerBound = cell(1);
			legend_names = cell(1);
			plot_num = 1;
		end

		for algo_i = 1:length(baseline_algos)
			algo = baseline_algos{algo_i}{1};
			algo_name = baseline_algos{algo_i}{2};
			known_algo = sum(ismember(baseline_algos_names, algo_name)) == 1;
			if (known_data_piece == false || known_algo == false)
				[curr_labels, curr_energy, curr_lowerBound, curr_time] = algo(unary, vertC, horC);
				labels{plot_num} = curr_labels;
				energy{plot_num} = curr_energy;
				lowerBound{plot_num} = curr_lowerBound;
				time{plot_num} = curr_time;
				step{plot_num} = [];
				oracle_calls{plot_num} = [];
				legend_names{plot_num} = algo_name;
				plot_num = plot_num + 1;
			end
		end
		% Set ylimits for making small size plots
		small_upper = energy{end}(5);
		small_lower = lowerBound{end}(5);

		% Now test all TRW variants
		for algo_i = 1:length(step_algos)
			algo = step_algos{algo_i}{1};
			init_params = step_algos{algo_i}{2};
			algo_name = step_algos{algo_i}{3};
			known_algo = sum(ismember(step_algos_names, algo_name)) == 1;
			if (known_data_piece == false || known_algo == false)
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
			if (~isempty(step{i}))
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
		figureProperties(plot_arr, curr_legend_names, 'Step plot', 'Time (sec)', 'Step size');
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
			if (~isempty(oracle_calls{i}))
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


		out_filename = strcat(experiment_folder, '/data_', data_name, '.mat');
		save(out_filename, 'labels', 'energy', 'oracle_calls', 'time', 'step', 'lowerBound', 'legend_names');
	end

	% Store what methods & data we used here.
	% This make possible to recognize already
	% processed data and not recompute it in future.
	data_set_mat = cell2mat(data_set);
	data_set_names = {data_set_mat.name};
	baseline_algos_mat = [baseline_algos{:}];
	baseline_algos_names = {baseline_algos_mat(2:2:end)};
	baseline_algos_names = baseline_algos_names{1};
	step_algos_mat = [step_algos{:}];
	step_algos_names = {step_algos_mat(3:3:end)};
	step_algos_names = step_algos_names{1};
	out_filename = strcat(experiment_folder, '/processed_list.mat');
	save(out_filename, 'data_set_names', 'baseline_algos_names', 'step_algos_names');
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
