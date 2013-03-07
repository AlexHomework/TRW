function compareEverything(experiment_name, varargin)
	% Compare different optimization methods, draw aggregate plots
	% and save everything (plots and the .mat file)
	% to experiments/experiment_name folder.
	% If you use this function with existing experiment folder,
	% it will try to use as much stored information as
	% possible to reduce computations.
	% 
	% Optional parameters:
	% 	profile_plots is a vector with iteration of TRW numbers
	%	 	on which profile plots (plot of the function
	%		in the optimization direction) will be shown.
	% 

	[profile_plots] = process_options(varargin, 'profile_plots', []);

	experiment_folder = strcat('experiments/', experiment_name);


	oracle_calls_counter = 0;
	function [value, der] = gridDualCounted(lambda, unary, vertC, horC)
		% Wrapper to count gridDual function calls
		oracle_calls_counter = oracle_calls_counter + 1;
		[value, der] = gridDual(lambda, unary, vertC, horC);
	end
	function counter = getOracleCalls()
		% Return oracle calls since last query or reset.
		% 
		counter = oracle_calls_counter;
		resetOracleCalls();
	end
	function counter = resetOracleCalls()
		oracle_calls_counter = 0;
	end


	data_set = getDataSet();
	baseline_algos = getBaselineAlgos();
	step_algos = getStepComputingAlgos();
	allAlgoCount = numel(baseline_algos) + numel(step_algos);
	colorsList = colorScheme(allAlgoCount);
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
		% is_known_data_piece == true iff current data piece
		% was already processed with some algorithms.
		is_known_data_piece = any(ismember(data_set_names, data_name));

		if (is_known_data_piece)
			results_filename = strcat(experiment_folder, '/data_', data_name, '.mat');
			load(results_filename, 'labels', 'energy', 'oracle_calls', 'time', 'step', ...
							'lowerBound', 'legend_names', 'unary', 'vertC', 'horC', ...
							'data_name', 'profile_info');
			plot_num = length(labels) + 1;
		else
			% Start everything from scratch.
			labels = {};
			energy = {};
			oracle_calls = {};
			time = {};
			step = {};
			lowerBound = {};
			profile_info = {};
			legend_names = {};
			plot_num = 1;
		end

		for algo_i = 1:length(baseline_algos)
			algo = baseline_algos{algo_i}{1};
			algo_name = baseline_algos{algo_i}{2};
			is_known_algo = any(ismember(baseline_algos_names, algo_name));
			if (is_known_data_piece == false || is_known_algo == false)
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
			is_known_algo = any(ismember(step_algos_names, algo_name));
			if (is_known_data_piece == false || is_known_algo == false)
				[curr_labels, curr_energy, curr_lowerBound, ...
							curr_time, curr_step, curr_oracle_calls, curr_profile_info] = ...
							trwGridPotts(unary, vertC, horC, algo, init_params, ...
											'save_iterations', profile_plots);
				labels{plot_num} = curr_labels;
				energy{plot_num} = curr_energy;
				lowerBound{plot_num} = curr_lowerBound;
				time{plot_num} = curr_time;
				step{plot_num} = curr_step;
				oracle_calls{plot_num} = curr_oracle_calls;
				profile_info{algo_i} = curr_profile_info;
				legend_names{plot_num} = step_algos{algo_i}{3};
				plot_num = plot_num + 1;
			end
		end

		plot_num = plot_num - 1;



		time_fig = plotAll(time, {energy, lowerBound}, legend_names, ...
									@(time, energy) ~isempty(energy), 'Time plot', ...
									'Time (sec)', 'Energy', colorsList, 'MultipleParts', true);
		out_filename = strcat('time_', data_name);
		saveFigure(time_fig, out_filename, experiment_folder);
		out_filename = strcat(out_filename, '_small');
		saveFigure(time_fig, out_filename, experiment_folder, 'ylim', [small_lower, small_upper]);
		close;



		outer_iterations_fig = plotAll([], {energy, lowerBound}, legend_names, ...
									@(time, energy) length(energy) > 1, 'Comparative plot', ...
									'Outer iterations', 'Energy', colorsList, 'MultipleParts', true);
		out_filename = strcat('outer_iterations_', data_name);
		saveFigure(outer_iterations_fig, out_filename, experiment_folder);
		out_filename = strcat(out_filename, '_small');
		saveFigure(outer_iterations_fig, out_filename, experiment_folder, 'ylim', [small_lower, small_upper]);
		close;


		step_fig = plotAll(time, step, legend_names, ...
									@(time, step) ~isempty(step), 'Step plot', ...
									'Time (sec)', 'Step size', colorsList);
		out_filename = strcat('step_', data_name);
		saveFigure(step_fig, out_filename, experiment_folder);
		close;


		oracle_fig = plotAll(oracle_calls, {energy, lowerBound}, legend_names, ...
									@(oracle_calls, energy) ~isempty(oracle_calls), 'Oracle calls plot', ...
									'Oracle calls', 'Energy', colorsList, 'MultipleParts', true);
		out_filename = strcat('oracle_', data_name);
		saveFigure(oracle_fig, out_filename, experiment_folder);
		out_filename = strcat(out_filename, '_small');
		saveFigure(oracle_fig, out_filename, experiment_folder, 'ylim', [small_lower, small_upper]);
		close;

		
		% Draw and save profile plots
		resolution = 100; % Points of profile plot count
		for main_algo_i = 1:length(step_algos)
			for piece = 1:length(profile_info{main_algo_i})
				curr_info = profile_info{main_algo_i}{piece};
				f = @(step) gridDualCounted(curr_info.lambda_first + step * curr_info.direction, ...
																		unary, vertC, horC);


				profile_fig = figure;
				hold all;

				resetOracleCalls();
				curr_legend_names = {};
				plot_arr = [];
				steps_arr = [];
				for algo_i = 1:length(step_algos)
					algo = step_algos{algo_i}{1};
					init_params = step_algos{algo_i}{2};
					algo_name = step_algos{algo_i}{3};
					if (strcmp(algo_name, step_algos{main_algo_i}{3}))
						curr_step = curr_info.step;
						f_val = f(curr_step);
						curr_oracle_calls = curr_info.oracle_calls;
					else
						[~, grad] = f(0);
						resetOracleCalls();
						[~, curr_step, f_val] = algo(f, curr_info.direction(:), grad(:), ...
															[curr_info.dual_energy], 1, init_params);
						curr_oracle_calls = getOracleCalls();
					end
					steps_arr(end + 1) = curr_step;
					p = plot(curr_step, f_val, '-o');
					plot_arr = [plot_arr, p];
					curr_legend_names{end + 1} = strcat(algo_name, ' (', ...
										int2str(curr_oracle_calls), ')');
				end

				% Plot from zero to two times the average step size
				end_point = 2 * mean(steps_arr);
				moments = [0:(end_point / resolution):end_point];
				values = arrayfun(f, moments);
				p = plot(moments, values);

				title = strcat(step_algos{main_algo_i}{3}, ', iteration ', int2str(curr_info.iteration));
				setFigureProperties(plot_arr, curr_legend_names, title, ...
											'Step size', 'Energy', 'legend_location', 'South');
				hold off;

				out_filename = strcat(experiment_folder, '/profile_', ...
											step_algos{main_algo_i}{3}, '_', ...
											int2str(curr_info.iteration), ...
											'_', data_name);
				set(gcf, 'PaperPositionMode', 'auto');
				print('-depsc2', strcat(out_filename, '.eps'));
				saveas(profile_fig, out_filename,'fig');
				close;
			end
		end


		out_filename = strcat(experiment_folder, '/data_', data_name, '.mat');
		save(out_filename, 'labels', 'energy', 'oracle_calls', 'time', 'step', ...
							'lowerBound', 'legend_names', 'unary', 'vertC', 'horC', ...
							'data_name', 'profile_info');

		% Store what methods & data we used here.
		% This make possible to recognize already
		% processed data and not recompute it in future.
		storeProcessedList(data_set(1:data_piece_ind), baseline_algos, step_algos, experiment_folder);
	end
end

function storeProcessedList(proc_data_set, proc_baseline_algos, proc_step_algos, experiment_folder)
	% Store 'that we've done' list.
	% 

	data_set_mat = cell2mat(proc_data_set(:));
	data_set_names = {data_set_mat.name};
	baseline_algos_mat = [proc_baseline_algos{:}];
	baseline_algos_names = {baseline_algos_mat(2:2:end)};
	baseline_algos_names = baseline_algos_names{1};
	step_algos_mat = [proc_step_algos{:}];
	step_algos_names = {step_algos_mat(3:3:end)};
	step_algos_names = step_algos_names{1};
	out_filename = strcat(experiment_folder, '/processed_list.mat');
	save(out_filename, 'data_set_names', 'baseline_algos_names', 'step_algos_names');
end

function saveFigure(figure_handler, name, experiment_folder, varargin)
	% Store figure in .eps and .fig formats
	% 
	% Optional parameters:
	% 	ylim
	% 

	[user_ylim] = process_options(varargin, 'ylim', []);
	if (~isempty(user_ylim))
		ylim(user_ylim);
	end

	out_filename = strcat(experiment_folder, '/', name);
	set(gcf, 'PaperPositionMode', 'auto');
	print('-depsc2', strcat(out_filename, '.eps'));
	saveas(figure_handler, out_filename,'fig');
end

function figure_handler = plotAll(x_data, y_data, legend_names, pickFunc, title_str, ...
																xlabel_str, ylabel_str, ...
																colorsList, varargin)
	% Draw one comparative plot.
	% It will draw only i-th plots such that pickFunc(x_data{i}, y_data{i}) == true
	% 
	% Optional parameters:
	% 	Extend [true] -- if true each plot will
	% 		be extended to the right with horizontal line
	% 	MultipleParts [false] -- if true y_data considered as cell array of plot parts
	% 	all setFigureProperties optional parameters
	% 

	[extend, multiple_plot_parts, unused_options] = process_options(varargin, ...
												'Extend', true, 'MultipleParts', false);

	if multiple_plot_parts
		plots_count = length(y_data{1});
	else
		plots_count = length(y_data);
	end


	if (extend)
		max_x = -Inf;
		for i = 1:length(x_data)
			x_curr = x_data{i};
			if multiple_plot_parts
				y_curr = y_data{1}{i};
			else
				y_curr = y_data{i};
			end
			if (pickFunc(x_curr, y_curr))
				if x_curr(end) > max_x
					max_x = x_curr(end);
				end
			end
		end
	end

	figure_handler = figure;
	hold all;
	plot_arr = [];
	curr_legend_names = {};
	for i = 1:plots_count
		if isempty(x_data)
			x_curr = [];
		else
			x_curr = x_data{i};
		end
		if multiple_plot_parts
			% There are multiple plot parts for each type of plot.
			y_curr = y_data{1}{i};
		else
			y_curr = y_data{i};
		end

		if (pickFunc(x_curr, y_curr))
			if (isempty(x_curr))
				x_curr = [1:length(y_curr)];
			end
			curr_color = colorsList(i, :);
			p = plot(x_curr, y_curr, 'Color', curr_color);
			setLineProperties(p);
			plot_arr(end + 1) = p;

			if (extend)
				x_extension = (x_curr(end) + 1) : max_x;
				y_extension = y_curr(end) * ones(size(x_extension));
				p = plot(x_extension, y_extension, '--', 'Color', curr_color);
				setLineProperties(p);
			end

			if multiple_plot_parts
				for plot_part = 2:length(y_data)
					y_curr = y_data{plot_part}{i};
					p = plot(x_curr, y_curr, 'Color', curr_color);
					setLineProperties(p);

					if (extend)
						y_extension = y_curr(end) * ones(size(x_extension));
						p = plot(x_extension, y_extension, '--', 'Color', curr_color);
						setLineProperties(p);
					end
				end
			end

			curr_legend_names{end + 1} = legend_names{i};
		end
	end
	setFigureProperties(plot_arr, curr_legend_names, title_str, xlabel_str, ylabel_str, unused_options{:});
	hold off;
end
