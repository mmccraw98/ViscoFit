classdef SSE_Fitting_Simultaneous
    %SSE_Fitting Class Containing All Info Necessary to Calculate the SSE
    %   This class, once initialized with the appropriate information, can
    %   calculate the Sum of Squared Errors (SSE) for multiple experiments
    %   at once.
    
    properties
        forces cell
        times cell
        indentations cell
        radii cell
        minTimescale double
    end
    
    methods
        function obj = SSE_Fitting_Simultaneous(forces,times,indentations,radii,minTimescale)
            %SSE_Fitting Construct an instance of the SSE_Fitting class
            %   This class is utilized to perform the parameter
            %   optimization steps. It helps collect the relevant data for
            %   multiple load-level fitting (i.e. multiple approach
            %   velocities at the same point). The properties (forces,
            %   times, indentations, and radii) must be supplied as
            %   cell arrays, where each entry corresponds to a different
            %   experiment. This is necessary for situations where the
            %   experiments are different lengths, to avoid needing to
            %   zero-pad the inputs and pass matrices to this class.
            if nargin >= 5
                obj.forces = forces;
                obj.times = times;
                obj.indentations = indentations;
                obj.radii = radii;
                obj.minTimescale = minTimescale;
            else
                error('Not enough data provided to the class definition. Please supply at least the force, time, indentation, and radii.');
            end
        end
        
        function sse = SSE_Maxwell(obj,params)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            
            % Calculate test forces
            test_forces = cellfun(@(t,h,r) LR_Maxwell(params,t,h,r),...
                obj.times, obj.indentations, obj.radii,...
                'UniformOutput',false);
            
            % calculate global residual
            sse_global = cellfun(@(F_data,F) sum((F_data-F).^2),...
                obj.forces, test_forces, 'UniformOutput',false);
            sse = sum(cell2mat(sse_global));
            
            if params(1) ~= 0
                elasticSetting = 'y';
            else
                elasticSetting = 'n';
            end
            if params(2) ~= 0
                fluidSetting = 'y';
            else
                fluidSetting = 'n';
            end
            
            [tauInds,modulusInds] = getParamIndices(params,elasticSetting,fluidSetting);
            ub = NaN(size(params));
            lb = NaN(size(params));
            
            ub(modulusInds) = 1e12;
            lb(modulusInds) = 1e-2;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end
        
        function sse = SSE_Voigt(obj,params)
            %SSE_Voigt Calculate the SSE for the Voigt model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Voigt Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            
            % Calculate test forces
            test_indentations = cellfun(@(t,F,r) LR_Voigt(params,t,F,r),...
                obj.times, obj.forces, obj.radii,...
                'UniformOutput',false);
            
            % calculate global residual
            sse_global = cellfun(@(h_data,h) sum(((h_data.^1.5)-h).^2),...
                obj.indentations, test_indentations, 'UniformOutput',false);
            sse = sum(cell2mat(sse_global));
            
            if params(1) ~= 0
                elasticSetting = 'y';
            else
                elasticSetting = 'n';
            end
            if params(2) ~= 0
                fluidSetting = 'y';
            else
                fluidSetting = 'n';
            end
            
            [tauInds,modulusInds] = getParamIndices(params,elasticSetting,fluidSetting);
            ub = NaN(size(params));
            lb = NaN(size(params));
            
            ub(modulusInds) = 1e2;
            lb(modulusInds) = 1e-12;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end
        
    end
end

