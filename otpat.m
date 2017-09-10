classdef otpat < handle
    properties (Access = public)
        cache
        measurement
        source
        load
        parameter
        model
        reg
    end
    
    methods
        function obj = otpat(opt)
            assert(isfield(opt, 'omega'));
            assert(isfield(opt, 'kappa'));
            assert(isfield(opt, 'femm_opt'));
            obj.cache = struct('e', [], 's', [], 'm', [], 'var', [], 'n', [], 'dof', [], 'ndof', [], 'k', 0.,  'o', 0.);
            obj.parameter = struct('d', [], 'a', [], 'g', []);
            obj.measurement = struct('h', [], 'j', []);
            obj.source = struct('pat', [], 'ot', []);
            obj.load = struct('pat', [], 'ot', []);
            
            obj.model = femm(opt.femm_opt);
            obj.cache.e = obj.model.build('e', 1, 'all');
            obj.cache.s = obj.model.build('s', 1);
            obj.cache.m = obj.model.build('m', 1);
            obj.cache.n = size(obj.model.space.nodes, 2);
            obj.cache.ndof = unique(obj.model.space.edges);
            obj.cache.dof  = setdiff(1:obj.cache.n, obj.cache.ndof);
            obj.cache.k = opt.kappa;
            obj.cache.o = opt.omega * 2 * pi/ (3e2);
            
            obj.reg.d = 1e0;
            obj.reg.a = 1e0;
            
            %% get PAT loads
            obj.source.pat = pat_source();
            obj.load.pat = zeros(obj.cache.n, length(obj.source.pat)); % allocation.
            
            for sId = 1:length(obj.source.pat)
                obj.load.pat(:, sId) = obj.model.build('g', obj.source.pat{sId}(obj.model.quad1d), 'all');
            end
            
            %% get OT loads
            obj.source.ot = ot_source();
            obj.load.ot = zeros(obj.cache.n, length(obj.source.ot)); % allocation.
            
            for sId = 1:length(obj.source.ot)
                obj.load.ot(:, sId) = obj.model.build('g', obj.source.ot{sId}(obj.model.quad1d), 'all');
            end
            
            %% get parameters 
            obj.parameter.a = absorptionF(obj.model.space.nodes)';
            obj.parameter.d = diffusionF(obj.model.space.nodes)';
            obj.parameter.g = gruneisenF(obj.model.space.nodes)';
            
        end
        
        
        function [m, U, V, A, B] = forward(obj, p, noise)
            if nargin == 2
                noise = 0;
            end
            %% get PAT, OT measurement
            qd = obj.mapping(p.d, obj.model.space.elems, obj.model.facet.ref');
            qa = obj.mapping(p.a, obj.model.space.elems, obj.model.facet.ref');

            S = obj.model.build('s', qd); 
            M = obj.model.build('m', qa);
            E = obj.cache.e / obj.cache.k;
            L = obj.load.pat / obj.cache.k;
            A =  (S + M + E); 
            U = A \ L;
            m.h = bsxfun(@times, (p.a .* p.g),U);
            m.h = m.h .* (1 + 2 * noise * (rand(size(m.h)) - 0.5));
            
            
            R = obj.load.ot / obj.cache.k;
            B = (S + M + E + sqrt(-1) * obj.cache.o * obj.cache.m);
            V = B \ R;
            m.j = V(obj.cache.ndof, :);
            m.j = m.j .* (1 + 2 * noise * (rand(size(m.j)) - 0.5));
            
        end
        
        
        function [f, g] = objective_gradient(obj, p)
            assert(size(p, 1) == 2 * obj.cache.n);
            local = struct('d', p(1:obj.cache.n),...
                'a', p(obj.cache.n+1:end), ...
                'g', ones(obj.cache.n, 1));
            [m, u, v, mp, mo] = forward(obj, local);
            z = obj.measurement;
            
            f = 0;        
            g_d = zeros(obj.cache.n, 1);
            g_a = zeros(obj.cache.n, 1);
            
            %% PAT part
            for sId = 1:size(m.h, 2)
                    
                crxRes = m.h(:, sId) - z.h(:,sId);
                f = f + 0.5 * obj.normsq(crxRes);

                g_a = g_a +  local.g .* (u(:, sId)) .* crxRes;

                %% adjoint
                phi =  mp \ (local.a .* local.g .* crxRes);

                g_d = g_d - ...
                    obj.model.adj(phi, u(:,sId), ones(obj.cache.n, 1), zeros(obj.cache.n, 1));
                    

                g_a = g_a - ...
                    obj.model.adj(phi, u(:,sId), zeros(obj.cache.n, 1), ones(obj.cache.n, 1));
            end
                    
            
            rd = obj.reg.d *(size(m.h, 2)) * obj.cache.n;
            ra = obj.reg.a *(size(m.h, 2)) * obj.cache.n;
            
            %% regularization
            f = f + 0.5 * rd *(local.d'* obj.cache.s * local.d) + 0.5 * ra * (local.a'*obj.cache.s*local.a);
            g_d = g_d + rd * obj.cache.s * local.d;
            g_a = g_a + ra * obj.cache.s * local.a;

            %% boundary prescribed
            g_d(obj.cache.ndof) = 0.;
            g_a(obj.cache.ndof) = 0.;

            g = [g_d; g_a];
             
        end
        
        function [res, hist] = backward(obj, init)
            if (nargin == 1)
                init = zeros(obj.cache.n * 2, 1);
                init(1:obj.cache.n) = 0.05;
                init(obj.cache.ndof) = obj.parameter.d(obj.cache.ndof);
                init(obj.cache.n + 1:end) = 0.1;
                init(obj.cache.n  + obj.cache.ndof) = obj.parameter.a(obj.cache.ndof);
                
            end
            opts    = struct( 'factr', 1e0, 'pgtol', 1e-10, 'm', 20, 'x0', init, 'maxIts', 600, 'maxTotalIts', 1e5);
            opts.printEvery     = 1;

            [res, ~, hist] =...
                lbfgsb_c(@obj.objective_gradient, zeros(size(init)), inf * ones(size(init)), opts); 

        end
        
        function vis(obj, h)
            trisurf(obj.model.space.elems(1:3,:)', obj.model.space.nodes(1,:), ...
            obj.model.space.nodes(2,:), h, 'EdgeColor','none');shading interp; colorbar;
            colormap jet;view(2);
        end
        
    end
    
    methods(Static)
        function [interpolate] = mapping(func, elems, trans_ref)
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(elems, 2));
            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end
        function r = normsq(v)
            r = sum(v.^2);
        end

    end
    
end

