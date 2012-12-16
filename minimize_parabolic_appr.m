function x_min = minimize_parabolic_appr(varargin)
	% Minimum point of parabolic approximation
	% Input examples:
	% 
	% 'd', x1, derivative(x1), 'v', x2, value(x2), v, x3, value(x3)
	% 
	% 'd', x1, derivative(x1), 'd', x2, derivative(x2)
	% 

	equations = zeros(3, 3);
	answers = zeros(3, 1);
	eq_num = 1;
	for i = 1:3:length(varargin)
		if strcmp(varargin{i}, 'd')
			x = varargin{i + 1};
			f_d_x = varargin{i + 2};
			equations(eq_num, :) = [2 * x, 1, 0];
			answers(eq_num) = f_d_x;
		elseif strcmp(varargin{i}, 'v')
			x = varargin{i + 1};
			f_x = varargin{i + 2};
			equations(eq_num, :) = [x^2, x, 1];
			answers(eq_num) = f_x;
		end
		eq_num = eq_num + 1;
	end
	if all(equations(3, :) == 0)
		% Only derivatives, special case
		x_min = minimize_der_parabolic_appr(varargin{2}, varargin{3}, varargin{5}, varargin{6});
	else
		coefficients = linsolve(equations, answers);
		a = coefficients(1);
		b = coefficients(2);
		x_min = -b / (2 * a);
	end
end
