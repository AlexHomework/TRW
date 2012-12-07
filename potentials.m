function [unary, vertC, horC] = potentials(name)
	filename = strcat('data/stereo/', name);
	left = RGB2Luv(imread(strcat(filename, '_l.ppm')));
	right = RGB2Luv(imread(strcat(filename, '_r.ppm')));
	[N, M, K] = size(left);

	vertC = 20 * ones(N - 1, M);
	horC = 20 * ones(N, M - 1);

	unary_filename = strcat(filename, '_datacost.txt');
	unary = readUnary(unary_filename, N, M);
end
