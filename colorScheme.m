function colorList = colorScheme(colorsCount)
	% Return vector of distinguishable colors to use for plotting.
	% 
	% Unfortunately, build in function
	% provide very similar colors (like 'hsv'),
	% or only few unique colors
	% available (like in the 'lines' function)
	% 

	colorList = [
		0.00  0.00  1.00;
		0.00  0.50  0.00;
		1.00  0.00  0.00;
		0.00  0.75  0.75;
		0.75  0.00  0.75;
		0.75  0.75  0.00;
		0.25  0.25  0.25;
		0.75  0.25  0.25;
		0.95  0.95  0.00;
		0.25  0.25  0.75;
		0.75  0.75  0.75;
		0.00  1.00  0.00;
		0.76  0.57  0.17;
		0.54  0.63  0.22;
		0.34  0.57  0.92;
		1.00  0.10  0.60;
		0.88  0.75  0.73;
		0.10  0.49  0.47;
		0.66  0.34  0.65;
		0.99  0.41  0.23;
	];

	if nargin > 0
		colorList = colorList(1:colorsCount, :);
	end
end