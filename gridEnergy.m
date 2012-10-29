function [energy] = gridEnergy(unary, vertC, horC, y)
	[K, N, M] = size(unary);
	% It's tricky vectorised equivalent to the code below
	energy = sum(unary(y(:) + K * (0 : N*M - 1)'));
	energy = energy + (y(1 : N*(M - 1)) ~= y(N+1 : end)) * horC(:);
	y_w_f = y(2:end, :); % y without first row
	y_w_l = y(1:N-1, :); % y without last row
	energy = energy + vertC(:)' * (y_w_l(:) ~= y_w_f(:));

	% energy = 0;
	% for i = 1:N
	% 	for j = 1:M
	% 		energy = energy + unary(y(i, j), i, j);
	% 		% Edges
	% 		if j > 1 && y(i, j - 1) ~= y(i, j)
	% 			energy = energy + horC(i, j - 1);
	% 		end
	% 		if i > 1 && y(i - 1, j) ~= y(i, j)
	% 			energy = energy + vertC(i - 1, j);
	% 		end
	% 	end
	% end
end
