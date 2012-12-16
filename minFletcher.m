function [alpha_min, f_min, status] = minFletcher(func, d, varargin)
% Default parameters
epsilon = 1e-5;
rho = 0.1;
sigma = 0.7;
tau = 0.1;
xi = 9;
max_iter = 100;
display = false;
for i = 1:2:length(varargin)
	if strcmp(varargin{i}, 'rho')
		rho = varargin{i + 1};
	elseif strcmp(varargin{i}, 'sigma')
		sigma = varargin{i + 1};
	elseif strcmp(varargin{i}, 'tau')
		tau = varargin{i + 1};
	elseif strcmp(varargin{i}, 'xi')
		xi = varargin{i + 1};
	elseif strcmp(varargin{i}, 'max_iter')
		max_iter = varargin{i + 1};
	elseif strcmp(varargin{i}, 'display')
		display = varargin{i + 1};
	end
end
		

num_oracle = 0;
function [value, derivative] = user_func(h)
	num_oracle = num_oracle + 1;
	[value, derivative] = func(h);
	derivative = derivative(:);
end

flag = 0;
a_l = 0;
a_u = 10^99;
[f_l, g] = user_func(0);
f_d_l = g' * d;
a_0 = 1;
for k = 1:max_iter
	first_ok = 0;
	if abs(a_l - a_0) < epsilon
		% Will be problems with parabolic approximations
		break;
	end
	[f_0, g] = user_func(a_0);
	if f_0 > f_l + rho * (a_0 - a_l) * f_d_l
		if a_0 < a_u
			a_u = a_0;
		end
		a_0_new = minimize_parabolic_appr('v', a_l, f_l, 'd', a_l, f_d_l, 'v', a_0, f_0);
		a_0_new = max(a_0_new, a_l + tau * (a_u - a_l));
		a_0_new = min(a_0_new, a_u - tau * (a_u - a_l));
		a_0 = a_0_new;
		[f_0, g] = user_func(a_0);
	else
		first_ok = 1;
	end
	f_d_0 = g' * d;
	if f_d_0 < sigma * f_d_l
		a_0_new = minimize_parabolic_appr('d', a_l, f_d_l, 'd', a_0, f_d_0);
		if a_0_new < a_0 + tau * (a_0 - a_l)
			a_0_new = a_0 + tau * (a_0 - a_l);
		elseif a_0_new > a_0 + xi * (a_0 - a_l)
			a_0_new = a_0 + xi * (a_0 - a_l);
		end
		a_l = a_0;
		a_0 = a_0_new;
		f_l = f_0;
		f_d_l = f_d_0;
	elseif first_ok
		% For current a_0 both both conditions are satisfied
		flag = 1;
		break;
	end
	% a_l- a_0
	% f_d_0 == f_d_l

	if display
		disp(['Iteration #', num2str(k), ...
			', current minimum: f(', num2str(a_0), ') = ', num2str(f_0)]);
	end
end
alpha_min = a_0;
f_min = f_0;
status = struct('flag', flag, 'num_oracle', num_oracle);
end
