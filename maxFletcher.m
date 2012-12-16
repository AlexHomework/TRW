function [alpha_max, f_max, status] = maxFletcher(func, d, varargin)
	function [value, derivative] = minus_func(h)
		[value, derivative] = func(h);
		value = -value;
		derivative = -derivative;
	end

	[alpha_min, f_min, status] = minFletcher(@minus_func, d, varargin);
	alpha_max = alpha_min;
	f_max = -f_min;
end