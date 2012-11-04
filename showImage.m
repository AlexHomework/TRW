function showImage(disparity)
	% figure;
	% a = 1;
	% depth = a ./ (disparity - 11) + a;
	% depth((depth > 2)) = 2;
	% min(min(depth))
	% imshow(mat2gray(depth));
	imshow(mat2gray(disparity));
end
