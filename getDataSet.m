function data_set = getDataSet()
	% Return cell array of potentials
	% data_set{i} = struct('unary', 'vertC', 'horC', 'name');
	% 

	file_names = {...
	'tsukuba', ...
	'venus',...
	'cones',...
	'art',...
	};


	data_set = {};
	for i = 1:length(file_names)
		[unary, vertC, horC] = potentials(file_names{i});
		data_set{i} = struct('unary', unary, 'vertC', vertC, 'horC', horC, 'name', file_names{i});
	end

end
