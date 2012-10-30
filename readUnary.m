function [unary] = read_unary(name, N, M)
	% Read data file in Alahari format
	filename = strcat('data/', name);
	
	file_id = fopen(filename, 'r');
	unary = fscanf(file_id, '%d');
	fclose(file_id);
	K = length(unary) / (N * M);
	unary_ = reshape(unary, M, N, K);

	unary = zeros(K, N, M);
	for i = 1:K
		unary(i, :, :) = reshape(unary_(:, :, i), M, N)';
	end
end
