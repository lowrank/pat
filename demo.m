N = 256;
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1));
nodes = [cos(theta); sin(theta)];

femm_opt = struct('deg', 1, 'qdeg', 2, 'min_area', 1e-3, 'edge', nodes);
opt = struct('femm_opt', femm_opt, 'omega', 2, 'kappa', 2.0);

o = otpat(opt);
[o.measurement, ~, ~, ~, ~] = o.forward(o.parameter);

[res, hist] = o.backward();



