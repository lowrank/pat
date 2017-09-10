N = 128;
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1));
nodes = [cos(theta); sin(theta)];

femm_opt = struct('deg', 2, 'qdeg', 6, 'min_area', 5e-4, 'edge', nodes);
opt = struct('femm_opt', femm_opt, 'omega', 2, 'kappa', 2.0);

obj = otpat(opt);
[obj.measurement, ~, ~, ~, ~] = obj.forward(obj.parameter, 0.01);

tic;
obj.reg = struct('d' , 1e-2, 'a', 1e-6);
[res, hist] = obj.backward();
toc;




