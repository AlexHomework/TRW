function [disparity] = stereo(name)
	filename = strcat('datasets/', name);
	left = RGB2Luv(imread(strcat(filename, '1.png')));
	right = RGB2Luv(imread(strcat(filename, '2.png')));
	[unary, vertC, horC] = potentials(left, right);
	labels = trwGridPotts(unary, vertC, horC);
end
