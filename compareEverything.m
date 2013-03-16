function compareEverything(experiment_name)
	% Compare different optimization methods, draw aggregate plots
	% and save everything (plots and the .mat file)
	% to experiments/experiment_name folder.
	% If you use this function with existing experiment folder,
	% it will try to use as much stored information as
	% possible to reduce computations.
	% 

	experiment_folder = strcat('experiments/', experiment_name);

	data_set = getDataSet();
	current_algos = getAlgos();
	if exist(experiment_folder, 'dir')
		% Find that we have already computed
		processed_filename = strcat(experiment_folder, '/data_sets.mat');
		load(processed_filename, 'data_set_names');
		processed_data_set =  data_set_names;
	else
		processed_data_set = {};
		mkdir(experiment_folder);
	end

	for data_piece_ind = 1:length(data_set)
		data_piece = data_set{data_piece_ind};
		unary = data_piece.unary;
		vertC = data_piece.vertC;
		horC = data_piece.horC;
		data_name = data_piece.name;
		% is_processed_data_piece == true iff current data piece
		% was already processed with some algorithms.
		is_processed_data_piece = any(ismember(processed_data_set, data_name));

		if (is_processed_data_piece)
			results_filename = strcat(experiment_folder, '/data_', data_name, '.mat');
			load(results_filename, 'labels', 'energy', 'oracle_calls', 'time', 'step', ...
							'lowerBound', 'algo_names', 'unary', 'vertC', 'horC', ...
							'data_name');
			plot_num = length(labels) + 1;
		else
			% Start everything from scratch.
			labels = {};
			energy = {};
			oracle_calls = {};
			time = {};
			step = {};
			lowerBound = {};
			algo_names = {};
			plot_num = 1;
		end

		for algo_i = 1:length(current_algos)
			curr_algo = current_algos{algo_i}{1};
			curr_algo_name = current_algos{algo_i}{2};
			is_processed_algo = any(ismember(algo_names, curr_algo_name));
			if (is_processed_data_piece == false || is_processed_algo == false)
				[curr_labels, curr_energy, curr_lowerBound, curr_time, curr_oracle_calls] = ...
																curr_algo(unary, vertC, horC);
				labels{plot_num} = curr_labels;
				energy{plot_num} = curr_energy;
				lowerBound{plot_num} = curr_lowerBound;
				time{plot_num} = curr_time;
				oracle_calls{plot_num} = curr_oracle_calls;
				algo_names{plot_num} = curr_algo_name;
				plot_num = plot_num + 1;
			end
		end
		plot_num = plot_num - 1;

		plot_algo_list = [];
		for algo_i = 1:length(algo_names)
			curr_algo_name = algo_names{algo_i};
			plot_algo_list(end + 1) = plotAlgo(current_algos, algo_names{algo_i});
		end
		% Set ylimits for making small size plots
		small_upper = energy{end}(5);
		small_lower = lowerBound{end}(5);


		colors_list = colorScheme(length(algo_names));
		time_fig = plotAll(time, {energy, lowerBound}, algo_names, plot_algo_list, ...
									@(time, energy) ~isempty(energy), 'Time plot', ...
									'Time (sec)', 'Energy', colors_list, 'MultipleParts', true);
		out_filename = strcat('time_', data_name);
		saveFigure(time_fig, out_filename, experiment_folder);
		out_filename = strcat(out_filename, '_small');
		saveFigure(time_fig, out_filename, experiment_folder, 'ylim', [small_lower, small_upper]);
		close;



		outer_iterations_fig = plotAll([], {energy, lowerBound}, algo_names, plot_algo_list, ...
									@(time, energy) length(energy) > 1, 'Comparative plot', ...
									'Outer iterations', 'Energy', colors_list, 'MultipleParts', true);
		out_filename = strcat('outer_iterations_', data_name);
		saveFigure(outer_iterations_fig, out_filename, experiment_folder);
		out_filename = strcat(out_filename, '_small');
		saveFigure(outer_iterations_fig, out_filename, experiment_folder, 'ylim', [small_lower, small_upper]);
		close;


		oracle_fig = plotAll(oracle_calls, {energy, lowerBound}, algo_names, plot_algo_list, ...
									@(oracle_calls, energy) ~isempty(oracle_calls), 'Oracle calls plot', ...
									'Oracle calls', 'Energy', colors_list, 'MultipleParts', true);
		out_filename = strcat('oracle_', data_name);
		saveFigure(oracle_fig, out_filename, experiment_folder);
		out_filename = strcat(out_filename, '_small');
		saveFigure(oracle_fig, out_filename, experiment_folder, 'ylim', [small_lower, small_upper]);
		close;


		out_filename = strcat(experiment_folder, '/data_', data_name, '.mat');
		save(out_filename, 'labels', 'energy', 'oracle_calls', 'time', 'step', ...
							'lowerBound', 'algo_names', 'unary', 'vertC', 'horC', ...
							'data_name');

		% Store list of processed data sets.
		out_filename = strcat(experiment_folder, '/data_sets.mat');
		data_set_mat = cell2mat(data_set(1:data_piece_ind));
		data_set_names = {data_set_mat.name};
		save(out_filename, 'data_set_names');
	end
end

function need_to_plot_algo = plotAlgo(algos_for_plotting, algo_name)
	% Check if we need to plot algo with name algo_name.
	% 

	need_to_plot_algo = false;
	for algo_i = 1:length(algos_for_plotting)
		if strcmp(algo_name, algos_for_plotting{algo_i}{2})
			need_to_plot_algo = true;
			break;
		end
	end
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

function figure_handler = plotAll(x_data, y_data, legend_names, plot_algo_list, pickFunc, ...
															title_str, xlabel_str, ylabel_str, ...
															colors_list, varargin)
	% Draw one comparative plot.
	% It will draw only i-th plots such that pickFunc(x_data{i}, y_data{i}) == true
	% and plot_algo_list(i) == true
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
		if plot_algo_list(i)
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
				curr_color = colors_list(i, :);
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
	end
	setFigureProperties(plot_arr, curr_legend_names, title_str, xlabel_str, ylabel_str, unused_options{:});
	hold off;
end
