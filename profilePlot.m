function fig = profilePlot(func, lambda_from, lambda_to, varargin)
	% Draw profile plot.
	% 
	% Optional arguments:
	% 	res [20] -- plot resolution
	%   title ['Profile plot']
	% 

	[resolution, title] = process_options(varargin, 'res', 20, 'title', 'Profile plot');
	from = 0;
	to = 2;
	moments = [from:((to - from) / resolution):to];
	direction = lambda_to - lambda_from;
	values = arrayfun(@(x) func(lambda_from + x * direction), moments);
	fig = figure;
	hold on;
	p1 = plot(moments, values, '-r');
	setLineProperties(p1);
	setFigureProperties(p1, {'Energy(base + step * direction)'}, title, 'Step', 'Energy', 'legend_location', 'South');
end