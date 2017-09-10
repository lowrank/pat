% unit disk
N = 128;
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1));
nodes = [cos(theta); sin(theta)];

femm_opt = struct('deg', 2, 'qdeg', 6, 'min_area', 2e-3, 'edge', nodes);
opt = struct('femm_opt', femm_opt, 'omega', 2, 'kappa', 2.0);

obj = otpat(opt);
[obj.measurement, ~, ~, ~, ~] = obj.forward(obj.parameter, 0.01);

tic;
obj.reg = struct('d' , 1e-2, 'a', 1e-6);
[res, hist] = obj.backward();
toc;
%%
m = 3; n = 3; order = 0;

set(gcf, 'Position', [100, 100, 1200, 800])
order = order + 1;
subplot(m, n, order);
obj.vis(obj.parameter.d);
title('exact  D');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(obj.parameter.a);
title('exact  \sigma');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(obj.parameter.d./(obj.parameter.a .* obj.parameter.g).^2);
title('exact \chi');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(res(1:obj.cache.n));
title('recovered D');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(res(obj.cache.n+1:2*obj.cache.n));
title('recovered \sigma');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(res(1:obj.cache.n)./(res(obj.cache.n+1:2*obj.cache.n)).^2);
title('recovered \chi');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(res(1:obj.cache.n) - obj.parameter.d);
title('error of D');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(res(obj.cache.n+1:2*obj.cache.n) - obj.parameter.a);
title('error of \sigma');
axis off;
order = order + 1;
subplot(m, n, order);
obj.vis(res(1:obj.cache.n)./(res(obj.cache.n+1:2*obj.cache.n)).^2 - obj.parameter.d./(obj.parameter.a.*obj.parameter.g).^2);
title('error of \chi');
axis off;




