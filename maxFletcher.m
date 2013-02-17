function [alpha_max, f_max, status] = maxFletcher(func, direction, f_0, g_0, varargin)
	function [value, derivative] = minus_func(h)
		[value, derivative] = func(h);
		value = -value;
		derivative = -derivative;
	end

	[alpha_min, f_min, status] = minFletcher(@minus_func, direction, f_0, g_0, varargin);
	alpha_max = alpha_min;
	f_max = -f_min;
end