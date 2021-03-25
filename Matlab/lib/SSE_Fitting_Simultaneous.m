classdef SSE_Fitting_Simultaneous
    %SSE_Fitting Class Containing All Info Necessary to Calculate the SSE
    %   This class, once initialized with the appropriate information, can
    %   calculate the Sum of Squared Errors (SSE) for multiple experiments
    %   at once.
    
    properties
        forces double
        times double {mustBePositive}
        dts double {mustBePositive}
        indentations double
        tipSize double {mustBePositive}
        tipGeom string
        minTimescale double {mustBePositive}
        nu double {mustBePositive}
    end
    
    methods
        function obj = SSE_Fitting_Simultaneous(forces,times,indentations,tipSize,minTimescale,varargin)
            %SSE_Fitting Construct an instance of the SSE_Fitting class
            %   This class is utilized to perform the parameter
            %   optimization steps. It helps collect the relevant data for
            %   multiple load-level fitting (i.e. multiple approach
            %   velocities at the same point). The properties (forces,
            %   times, indentations, and tip size) must be supplied as
            %   cell arrays, where each entry corresponds to a different
            %   experiment. This is necessary for situations where the
            %   experiments are different lengths, to avoid needing to
            %   zero-pad the inputs and pass matrices to this class.
            
            %   Note, that unless specified in the varargin, the tipSize is
            %   assumed to be a tip radius (spherical indentation). If the
            %   second varargin is set to 'conical', then the contact
            %   mechanics are changed to account for the tipSize being the
            %   tip angle.
            
            if nargin >= 5
                obj.forces = horzcat(forces{:});
                obj.times = horzcat(times{:});
                obj.indentations = horzcat(indentations{:});
                temp = cellfun(@(t) mode(gradient(t)).*ones(size(t)),times,'UniformOutput',false);
                obj.dts = horzcat(temp{:});
                temp = cellfun(@(r,t) r.*ones(size(t)),tipSize,times,'UniformOutput',false);
                obj.tipSize = horzcat(temp{:});
                obj.minTimescale = minTimescale;
            else
                error('Not enough data provided to the class definition. Please supply at least the force, time, indentation, and radii.');
            end
            
            % Handle optional arguments. Create default values, and modify
            % if the varargin contains the designated "real" values.
            obj.nu = cell(size(forces));
            obj.nu(:) = 0.5;
            temp = cellfun(@(nu,t) nu.*ones(size(t)),obj.nu,times,'UniformOutput',false);
            obj.nu = horzcat(temp{:});
            obj.tipGeom = "spherical";
            for i=1:length(varargin)
                switch(i)
                    case 1
                        temp = cellfun(@(nu,t) nu.*ones(size(t)),varargin{i},times,'UniformOutput',false);
                        obj.nu = horzcat(temp{:});
                    case 2
                        if ~strcmp(obj.tipGeom, string(varargin{i}))
                            obj.tipGeom = string(varargin{i});
                        end
                end
            end
        end
        
        function sse = SSE_Maxwell(obj,params)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            
            % Calculate test forces
            test_forces = LR_Maxwell(params,obj.times,obj.dts,obj.indentations,obj.tipSize,obj.nu,obj.tipGeom);
            
            % calculate global residual
            sse_global = sum((obj.forces-test_forces).^2);
            sse = sum(sse_global);
            
            [tauInds,modulusInds] = getParamIndices(params);
            ub = zeros(size(params))+eps;
            lb = zeros(size(params));
            
            ub(modulusInds) = 1e12;
            lb(modulusInds) = 1e-2;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;
            
            if params(2) > 2*eps
                ub(2) = max(tauCenters)*1e2;
                lb(2) = min(obj.dts);
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end
        
        function sse = SSE_Voigt(obj,params)
            %SSE_Voigt Calculate the SSE for the Voigt model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Voigt Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            
            % Calculate test indentation
            test_indentations = LR_Voigt(params,obj.times,obj.dts,obj.forces,obj.tipSize,obj.nu,obj.tipGeom);
            
            switch obj.tipGeom
                case "spherical"
                    beta = 1.5;
                case "conical"
                    beta = 2;
            end
            
            % calculate global residual
            sse_global = sum(((obj.indentation.^beta)-test_indentations).^2);
            sse = sum(sse_global);
            
            [tauInds,modulusInds] = getParamIndices(params);
            ub = zeros(size(params))+eps;
            lb = zeros(size(params));
            
            ub(modulusInds) = 1e2;
            lb(modulusInds) = 1e-12;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;
            
            if params(2) > 2*eps
                ub(2) = 1;
                lb(2) = 0;
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end
        
    end
end

