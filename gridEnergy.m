function [energy] = gridEnergy(unary, vertC, horC, metric, y)
	energy = 0;
	[K, N, M] = size(unary);
	for i = 1:N
		for j = 1:M
			energy = energy + unary(y(i, j), i, j);

			% Edges
			if j > 1
				energy = energy + horC(i, j - 1) * metric(y(i, j - 1), y(i, j));
			end
			if j < M
				energy = energy + horC(i, j) * metric(y(i, j), y(i, j + 1));
			end
			if i > 1
				energy = energy + vertC(i - 1, j) * metric(y(i - 1, j), y(i, j));
			end
			if i < N
				energy = energy + vertC(i, j) * metric(y(i, j), y(i + 1, j));
			end
		end
	end
end
