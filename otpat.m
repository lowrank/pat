classdef otpat < handle
    properties (Access = public)
        cache
        measurement
        source
        load
        parameter
        model
    end
    
    methods
        function obj = otpat(opt)
            assert(isfield(opt, 'omega'));
            assert(isfield(opt, 'kappa'));
            assert(isfield(opt, 'femm_opt'));
            obj.cache = struct('e', [], 's', [], 'm', [], 'var', [], 'n', []);
            obj.parameter = struct('d', [], 'a', [], 'g', []);
            obj.measurement = struct('h', [], 'j', []);
            obj.source = struct('pat', [], 'ot', []);
            obj.load = struct('pat', [], 'ot', []);
            
            obj.model = femm(opt.femm_opt);
            obj.cache.e = obj.model.build('e', 1, 'all');
            obj.cache.s = obj.model.build('s', 1);
            obj.cache.m = obj.model.build('m', 1);
            obj.cache.n = size(obj.model.space.nodes, 2);
            
            %% get PAT loads
            obj.source.pat = pat_source();
            obj.load.pat = zeros(obj.cache.n, length(obj.source.pat)); % allocation.
            
            for sId = 1:length(obj.source.pat)
                obj.load.pat(:, sId) = obj.model.build('g', obj.source.pat{sId}(obj.model.quad1d), 'all');
            end
            
            %% get OT loads
            
            
            
        end
        
        
        function forward(obj)
            %% get PAT, OT measurement
            
        end
        
        
        
        
    end
    
end

