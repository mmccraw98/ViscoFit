classdef ViscoFit
    %ViscoFit Class Containing All Info Necessary to Extract Viscoelastic
    %Parameters
    %   This class, once initialized with the appropriate information, can
    %   perform a variety of viscoelastic parameterization operations. It
    %   must be initialized with the observed force, time, indentation, tip
    %   size (either spherical tip radius or indenter cone angle), the
    %   minimum timescale to begin fitting with: 
    %   "forces" - {1xN} - Cell array containing repulsive force array from
    %   experiments 1 to N (if more than one experiment is considered). The
    %   units are NEWTONS. The arrays must use DOUBLE PRECISION.
    
    %   "times" - {1xN} - Cell array containing repulsive time array from
    %   experiments 1 to N (if more than one experiment is considered). The
    %   units are SECONDS. The arrays must use DOUBLE PRECISION.
    
    %   "indentations" - {1xN} - Cell array containing repulsive
    %   indentation observed for experiments 1 to N (if more than one 
    %   experiment is considered). The units are METERS. The arrays must 
    %   use DOUBLE PRECISION.
    
    %   "tipSize" - {1xN} - Cell array containing the characteristic tip
    %   size, either the tip radius (for spherical indentation, the
    %   default) or cone angle (for conical indentation). Note that for
    %   conical indentation, an additional argument must be passed to
    %   overwrite the default spherical indentation setting. The units are
    %   METERS (spherical) or DEGREES (conical).  The values must use 
    %   DOUBLE PRECISION.
    
    %   "minTimescale" - double - Single value which will be the "center"
    %   of the characteristic time range for the first viscoelastic element
    %   in the generalized rheological models. This is not necessary when
    %   using the PLR model; when using PLR, set this value to ONE.  The 
    %   value must use DOUBLE PRECISION.
    
    %   Additionally, the class can take a variety of other settings:
    
    %   Optional Argument #1 - the sample's poisson's ratio, if known, and
    %   not equal to 0.5 (incompressible) which is the default.
    
    %   Optional Argument #2 - the tip geometry, a string of either
    %   "spherical" or "conical", which determines how the tipSize argument
    %   is interpreted.
    
    properties
        forces double
        forces_cell cell
        times double {mustBePositive}
        times_cell cell
        dts double {mustBePositive}
        dts_cell cell
        indentations double
        indentations_cell cell
        tipSize double {mustBePositive}
        tipSize_cell cell
        tipGeom string
        minTimescale double {mustBePositive}
        nu double {mustBePositive}
        nu_cell cell
        fitLog logical
    end
    
    methods
        function obj = ViscoFit(forces,times,indentations,tipSize,varargin)
            %ViscoFit Construct an instance of the ViscoFit class
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

            if nargin >= 4
                
                % Handle optional arguments. Create default values, and modify
                % if the varargin contains the designated "real" values.
                tempnu = cell(size(forces));
                tempnu(:) = {0.5}; % Incompressible, 0.5
                temp = cellfun(@(nu,t) nu.*ones(size(t)),tempnu,times,'UniformOutput',false);
                obj.nu = horzcat(temp{:});
                obj.nu_cell = temp;
                obj.tipGeom = "spherical";
                obj.fitLog = false;
                obj.minTimescale = 1e-4; % Default
                if ~isempty(varargin)
                    if isa(varargin{1},'struct')
                        % Grab the settings structure
                        inputSettings = varargin{1};
                        
                        % Store the settings in the class
                        if ~logical(inputSettings.fitLog)
                            temp = cellfun(@(nu,t) nu.*ones(size(t)),inputSettings.nu,times,'UniformOutput',false);
                            obj.nu = horzcat(temp{:});
                            obj.nu_cell = temp;
                        else
                            tempdt = cellfun(@(t) round(mode(gradient(t)),1,'significant'),times,'UniformOutput',false);
                            temp = cellfun(@(x,t,dt)log_scale(x.*ones(size(t)),t,mode(dt),t(end)),inputSettings.nu,times,tempdt,'UniformOutput',false);
                            obj.nu = horzcat(temp{:});
                            obj.nu_cell = temp;
                        end
                        obj.tipGeom = string(inputSettings.tipGeom);
                        obj.fitLog = inputSettings.fitLog;
                        obj.minTimescale = inputSettings.minTimescale;
                        
                    else
                        error('You are not passing the settings correctly to the ViscoFit Class Initialization. Please ensure the fifth argument is a structure containing your settings.');
                    end
                    
                end
                
                if ~obj.fitLog
                    % Use the full-fidelity data (no log sampling)
                    obj.forces = horzcat(forces{:});
                    obj.forces_cell = forces;

                    obj.times = horzcat(times{:});
                    obj.times_cell = times;

                    obj.indentations = horzcat(indentations{:});
                    obj.indentations_cell = indentations;

                    temp = cellfun(@(t) round(mode(gradient(t)),1,'significant').*ones(size(t)),times,'UniformOutput',false);
                    obj.dts = horzcat(temp{:});
                    obj.dts_cell = temp;

                    temp = cellfun(@(r,t) r.*ones(size(t)),tipSize,times,'UniformOutput',false);
                    obj.tipSize = horzcat(temp{:});
                    obj.tipSize_cell = temp;
                else
                    % Log-sample all of the data before performing any
                    % operations. Note: this reduces the data quality, but
                    % drastically improves the fitting speed for legacy
                    % methods (NLS, in particular).
                    tempdt = cellfun(@(t) round(mode(gradient(t)),1,'significant'),times,'UniformOutput',false);
                    temp = cellfun(@(x,t,dt)log_scale(x,t,mode(dt),t(end)),times,times,tempdt,'UniformOutput',false);
                    obj.times = horzcat(temp{:});
                    obj.times_cell = temp;
                    
                    temp = cellfun(@(t,x) round(mode(gradient(t)),1,'significant').*ones(size(x)),times,obj.times_cell,'UniformOutput',false);
                    obj.dts = horzcat(temp{:});
                    obj.dts_cell = temp;
                    
                    temp = cellfun(@(x,t,dt)log_scale(x,t,mode(dt),t(end)),forces,times,obj.dts_cell,'UniformOutput',false);
                    obj.forces = horzcat(temp{:});
                    obj.forces_cell = temp;

                    temp = cellfun(@(x,t,dt)log_scale(x,t,mode(dt),t(end)),indentations,times,obj.dts_cell,'UniformOutput',false);
                    obj.indentations = horzcat(temp{:});
                    obj.indentations_cell = temp;
                    
                    temp = cellfun(@(r,t) r.*ones(size(t)),tipSize,obj.times_cell,'UniformOutput',false);
                    obj.tipSize = horzcat(temp{:});
                    obj.tipSize_cell = temp;
                    
                end
                
            else
                error('Not enough data provided to the class definition. Please supply at least the force, time, indentation, and radii.');
            end
            
        end
        
        function sse = SSE_Maxwell(obj,params,elasticSetting,fluidSetting)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting simultaneously for all force
            %   curves by linearizing the entire dataset (i.e. taking all
            %   curves for a particular approach velocity and stacking them
            %   end-to-end in a row vector). Separation and convolution of
            %   the correct regions of this row vector is handled inside
            %   the LR_Maxwell function.
            
            % Calculate test forces
            test_forces = LR_Maxwell(params,obj.times,obj.dts,obj.indentations,obj.tipSize,obj.nu,obj.tipGeom,elasticSetting,fluidSetting);
            
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
            
            if length(params) > 1
                if fluidSetting
                    ub(2) = max(tauCenters)*1e2;
                    lb(2) = min(obj.dts);
                end
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end % End Maxwell SSE Function
        
        function sse = SSE_Maxwell_Map(obj,params,idx,elasticSetting,fluidSetting)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting for all force curves
            %   separately, by taking in an additional index (idx) compared
            %   to the standard SSE_Maxwell function. This is intended
            %   solely for Force Map analysis, wherein each pixel
            %   (containing a single force curve) is treated for analysis.
            %   The output is the same as for the SSE_Maxwell function - a
            %   Sum of Squared Errors for that particular pixel.
            
            % Calculate test forces
            test_forces = LR_Maxwell(params,obj.times_cell{idx},obj.dts_cell{idx},obj.indentations_cell{idx},obj.tipSize_cell{idx},obj.nu_cell{idx},obj.tipGeom,elasticSetting,fluidSetting);
            
            % calculate global residual
            sse_global = sum((obj.forces_cell{idx}-test_forces).^2);
            sse = sum(sse_global);
            
            [tauInds,modulusInds] = getParamIndices(params);
            ub = zeros(size(params))+eps;
            lb = zeros(size(params));
            
            ub(modulusInds) = 1e12;
            lb(modulusInds) = 1e-2;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;
            
            if length(params) > 1
                if fluidSetting
                    ub(2) = max(tauCenters)*1e2;
                    lb(2) = min(obj.dts);
                end
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end % End Maxwell SSE Function
        
        function sse = SSE_Voigt(obj,params,elasticSetting,fluidSetting)
            %SSE_Voigt Calculate the SSE for the Voigt model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Voigt Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting simultaneously for all force
            %   curves by linearizing the entire dataset (i.e. taking all
            %   curves for a particular approach velocity and stacking them
            %   end-to-end in a row vector). Separation and convolution of
            %   the correct regions of this row vector is handled inside
            %   the LR_Voigt function.
            
            % Calculate test indentation
            test_indentations = LR_Voigt(params,obj.times,obj.dts,obj.forces,obj.tipSize,obj.nu,obj.tipGeom,elasticSetting,fluidSetting);
            
            switch obj.tipGeom
                case "spherical"
                    beta = 1.5;
                case "conical"
                    beta = 2;
            end
            
            % calculate global residual
            sse_global = sum(((obj.indentations.^beta)-test_indentations).^2);
            sse = sum(sse_global);
            
            [tauInds,modulusInds] = getParamIndices(params);
            ub = zeros(size(params))+eps;
            lb = zeros(size(params));
            
            ub(modulusInds) = 1e2;
            lb(modulusInds) = 1e-12;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;

            if length(params) > 1
                if fluidSetting
                    ub(2) = 1;
                    lb(2) = 0;
                end
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end % End Voigt SSE Map Function
        
        function sse = SSE_Voigt_Map(obj,params,idx,elasticSetting,fluidSetting)
            %SSE_Voigt Calculate the SSE for the Voigt model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Voigt Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting for all force curves
            %   separately, by taking in an additional index (idx) compared
            %   to the standard SSE_Voigt function. This is intended
            %   solely for Force Map analysis, wherein each pixel
            %   (containing a single force curve) is treated for analysis.
            %   The output is the same as for the SSE_Voigt function - a
            %   Sum of Squared Errors for that particular pixel.
            
            % Calculate test indentation
            test_indentations = LR_Voigt(params,obj.times_cell{idx},obj.dts_cell{idx},obj.forces_cell{idx},obj.tipSize_cell{idx},obj.nu_cell{idx},obj.tipGeom,elasticSetting,fluidSetting);
            
            switch obj.tipGeom
                case "spherical"
                    beta = 1.5;
                case "conical"
                    beta = 2;
            end
            
            % calculate global residual
            sse_global = sum(((obj.indentations_cell{idx}.^beta)-test_indentations).^2);
            sse = sum(sse_global);
            
            [tauInds,modulusInds] = getParamIndices(params);
            ub = zeros(size(params))+eps;
            lb = zeros(size(params));
            
            ub(modulusInds) = 1e2;
            lb(modulusInds) = 1e-12;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;

            if length(params) > 1
                if fluidSetting
                    ub(2) = 1;
                    lb(2) = 0;
                end
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end % End Voigt SSE Map Function
        
        function sse = SSE_PLR(obj,params,varargin)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting simultaneously for all force
            %   curves by linearizing the entire dataset (i.e. taking all
            %   curves for a particular approach velocity and stacking them
            %   end-to-end in a row vector). Separation and convolution of
            %   the correct regions of this row vector is handled inside
            %   the LR_PLR function.
            
            % Calculate test forces
            test_forces = LR_PLR(params,obj.times,obj.dts,obj.indentations,obj.tipSize,obj.nu,obj.tipGeom);
            
            % calculate global residual
            sse_global = sum((obj.forces-test_forces).^2);
            sse = sum(sse_global);
            
            % Power Law Rheology Roster:
            % [E_0 alpha]
            ub = [1e12;1];
            lb = [1e-2;0];
            
            if any(ub(1:length(params))-params < 0) || any(params-lb(1:length(params)) < 0)
                sse = Inf;
            end
        end % End PLR SSE Function
        
        function sse = SSE_PLR_Map(obj,params,idx,varargin)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting for all force curves
            %   separately, by taking in an additional index (idx) compared
            %   to the standard SSE_PLR function above. This is intended
            %   solely for Force Map analysis, wherein each pixel
            %   (containing a single force curve) is treated for analysis.
            %   The output is the same as for the SSE_PLR function - a
            %   Sum of Squared Errors for that particular pixel.
            
            % Calculate test forces
            test_forces = LR_PLR(params,obj.times_cell{idx},obj.dts_cell{idx},obj.indentations_cell{idx},obj.tipSize_cell{idx},obj.nu_cell{idx},obj.tipGeom);
            
            % calculate global residual
            sse_global = sum((obj.forces_cell{idx}-test_forces).^2);
            sse = sum(sse_global);
            
            % Power Law Rheology Roster:
            % [E_0 alpha]
            ub = [1e12;1];
            lb = [1e-2;0];
            
            if any(ub(1:length(params))-params < 0) || any(params-lb(1:length(params)) < 0)
                sse = Inf;
            end
        end % End PLR SSE Map Function
        
        function fitStruct = fitData(obj, varargin)
            % FITDATA Fit a Viscoelastic Model to the Class Data 
            %   This function takes in a variety of settings in addition to
            %   the data already provided to the class and performs an
            %   optimization procedure based on those settings. Defaults
            %   are set below, and the order in which they are provided
            %   must be adhered to (unless the user manually reprograms
            %   this class). In particular, the number of iterations
            %   (n_iterations) may need to be increased, as with the number
            %   of solver iterations (n_fitIterations) and maximum number
            %   of elements in the viscoelastic series (n_elements). This
            %   function will iteratively introduce viscoelastic elements,
            %   feeding forward the results from the previous iteration,
            %   until the optimization has been performed for a model
            %   configuration with a number of viscoelastic elements 
            %   equal to n_elements.
            
            % Initialize Output Structure
            fitStruct = struct;
            
            % Default Settings
            solver = 'nelder-mead';     % Fit using Nelder-Mead Simplex
            model = 'maxwell';          % Use Generalized Maxwell Model
            n_elements = 3;             % Fit iteratively for up to 3 elements
            elasticSetting = 1;         % Include Elastic Term
            fluidSetting = 0;           % No Steady-State Fluidity
            n_iterations = 100;         % Use 100 random initializations as a default
            n_fitIterations = 1e4;      % No. of iterations for solver
            if ~isempty(varargin)
                % Only one varargin is accepted, and it is a structure
                % containing all of the settings information we require
                fitOpts = varargin{1};
                
                % Get the settings provided
                solver = lower(fitOpts.solver);
                model = lower(fitOpts.model);
                n_elements = fitOpts.n_elements;
                elasticSetting = fitOpts.elasticSetting;
                fluidSetting = fitOpts.fluidSetting;
                n_iterations = fitOpts.n_iterations;
                n_fitIterations = fitOpts.n_fitIterations;
                if strcmpi(solver,'custom')
                    if isfield(fitOpts,'customFunc')
                        customFunc = fitOpts.customFunc;
                    else
                        error('You selected "custom" for the solver type and did not supply the function name in the class settings.');
                    end
                end
            end
            
            % Store the fit settings for future reference
            fitStruct.solver = solver;
            fitStruct.model = model;
            fitStruct.n_elements = n_elements;
            fitStruct.elasticSetting = elasticSetting;
            fitStruct.fluidSetting = fluidSetting;
            fitStruct.n_iterations = n_iterations;
            
            % Get the correct objective function for optimization
            switch model
                case 'maxwell'
                    objFunc = @obj.SSE_Maxwell;
                    modelFile = 'LR_Maxwell.m';
                case 'voigt'
                    objFunc = @obj.SSE_Voigt;
                    modelFile = 'LR_Voigt.m';
                case 'plr'
                    objFunc = @obj.SSE_PLR;
                    modelFile = 'LR_PLR.m';
                    % The storage roster for plr is different, and requires
                    % the second index to maintain consistency. Thus, we
                    % have to manually force the second term to be fit by
                    % including the fluidity. This second position actually
                    % corresponds to the exponent, alpha.
                    fluidSetting = 1;
                    fitStruct.fluidSetting = fluidSetting;
                case 'custom'
                    objFunc = @obj.customFunc;
                    modelFile = 'customFunc.m';
                otherwise
                    error('Your chosen solver-model combination is not implemented yet.');
            end
            
            % Create placeholders for the data we will obtain from our
            % optimization attempts
            fitStruct.bestParams = {};
            fitStruct.paramPopulation = {};
            fitStruct.paramPopulationResiduals = {};
            fitStruct.elasticFitTime = {};
            fitStruct.fitTime = {};
            
            % Open Parallel Pool of MATLAB Workers
            if isempty(gcp('nocreate'))
                % Make a fresh pool
                poolobj = parpool('IdleTimeout', Inf);
            else
                % Remove the old pool
                poolobj = gcp('nocreate');
                delete(poolobj);
                
                % Make a fresh pool
                poolobj = parpool('IdleTimeout', Inf);
            end
            
            % Send the class to the workers
            addAttachedFiles(poolobj, {'ViscoFit.m',modelFile})

            % Start the timer
            tic;
            
            % Begin the iterative term introduction loop
            for i = 1:n_elements   
                
                % Look to see if there are old results available to provide
                % intelligent guesses for our initial conditions. This will
                % only occur for iterations beyond the first
                if i > 1
                    % Our current loop has an array of parameters
                    % two-larger than the previous (for the generalized
                    % spring-dashpot rheology models)
                    beta_in = NaN(length(fitStruct.bestParams{i-1})+2,1);
                    
                    % Overwrite the NaN values with the previous optimal
                    % parameters.
                    beta_in(1:length(fitStruct.bestParams{i-1})) = fitStruct.bestParams{i-1};
                else
                    % For the first iteration, there will be four
                    % parameters in total: the elastic element (1), the
                    % steady-state fluidity (2), and two parameters for the
                    % first viscoelastic element (3,4). Note, based on the
                    % elasticSetting and fluidSetting provided to this
                    % function, the values in (1) and (2) may never be
                    % updated
                    if strcmpi(model,'plr')
                        beta_in = NaN(2,1);
                    else
                        beta_in = NaN(4,1);
                    end
                end
                
                % Make our storage variables
                beta_dist = zeros(length(beta_in),n_iterations);
                beta_dist_elastic = zeros(1,n_iterations);
                beta0_dist = zeros(size(beta_dist));
                residual_dist = NaN(1,n_iterations);
                residual_dist_elastic = NaN(size(residual_dist));
                
                % Make the upper and lower bounds for this model
                [tauInds,modulusInds] = getParamIndices(beta_in);
                ub = zeros(size(beta_in))+eps;
                lb = zeros(size(beta_in));
                switch model
                    case 'maxwell'
                        % Moduli are limited to "reasonable" bounds for
                        % viscoelastic materials
                        ub(modulusInds) = 1e12;
                        lb(modulusInds) = 1e-2;

                        tauCenters = obj.minTimescale.*(10.^( (1:length(ub(3:2:end)))-1 ));
                        ub(tauInds) = tauCenters*10;
                        lb(tauInds) = tauCenters/10;
                        
                        if ~elasticSetting
                            ub(1) = eps;
                            lb(1) = 0;
                        end

                        if fluidSetting
                            ub(2) = max(tauCenters)*1e2;
                            lb(2) = 10^( floor(min(log10(obj.dts)))+1 );
                        end
                        
                        % Restrict the range of random guesses, if desired.
                        % Otherwise, they should be set equal to ub & lb
                        ub_rand(modulusInds) = 6;
                        ub_rand(tauInds) = log10(ub(tauInds));
                        lb_rand(modulusInds) = 1;
                        lb_rand(tauInds) = log10(lb(tauInds));
                        
                    case 'voigt'
                        % Compliances are limited to "reasonable" bounds
                        % for viscoelastic materials
                        ub(modulusInds) = 1e2;
                        lb(modulusInds) = 1e-12;

                        tauCenters = obj.minTimescale.*(10.^( (1:length(ub(3:2:end)))-1 ));
                        ub(tauInds) = tauCenters*10;
                        lb(tauInds) = tauCenters/10;
                        
                        if ~elasticSetting
                            ub(1) = eps;
                            lb(1) = 0;
                        end

                        if fluidSetting
                            ub(2) = 1;
                            lb(2) = 0;
                        end

                        % Restrict the range of random guesses, if desired.
                        % Otherwise, they should be set equal to ub & lb
                        ub_rand(modulusInds) = -1;
                        ub_rand(tauInds) = log10(ub(tauInds));
                        lb_rand(modulusInds) = -6;
                        lb_rand(tauInds) = log10(lb(tauInds));
                        
                    case 'plr'
                        % Power Law Rheology Roster:
                        % [E_0 alpha]
                        ub = [1e12 1];
                        lb = [1e-2 0];
                        
                        % Restrict the range of random guesses, if desired.
                        % Otherwise, they should be set equal to ub & lb
                        ub_rand = log10(ub);
                        lb_rand = [log10(lb(1)) -3];
                        
                    case 'custom'
                        
                        % Define the upper and lower bounds for your custom
                        % function here. The output from this region should
                        % be two arrays, ub and lb, which are the same
                        % length as the number of parameters in the model.
                        % ub = [...];
                        % lb = [...];
                        
                end
                
                preFitting = toc;
                switch solver
                    case 'nelder-mead'
                        
                        options = optimset('Display','none',...
                                    'PlotFcns',[],...
                                    'MaxFunEvals',n_fitIterations,...
                                    'MaxIter',n_fitIterations,...
                                    'TolFun',1e-60,...
                                    'TolX',1e-60);
                        
                        if i == 1 && elasticSetting
                            % Fit the elastic term separately for the first
                            % iteration. Future iterations have the "best
                            % fit" elastic term included from the prior
                            % optimization attempt
                            
                            % Clock the timer
                            preElasticTime = toc;
                            
                            progressString = sprintf('Nelder-Mead (%s)\nInvestigating Elastic Parameter\nParallel Search Running...',model);
                            hbar = parfor_progressbar(n_iterations,progressString);
                            warning('off');
                            parfor j = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fminsearch(@(x)objFunc(x,elasticSetting,fluidSetting),beta0,options);
                                hbar.iterate(1);
                            end
                            close(hbar);
                            warning('on');
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = horzcat(fitStruct.elasticFitTime,{postElasticTime-preElasticTime});
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        progressString = sprintf('Nelder-Mead (%s)\nInvestigating Material Parameters\nParallel Search Running...',model);
                        hbar = parfor_progressbar(n_iterations,progressString);
                        warning('off');
                        parfor j = 1:n_iterations
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = fminsearch(@(x)objFunc(x,elasticSetting,fluidSetting),beta0,options);
                            hbar.iterate(1);
                        end
                        close(hbar);
                        warning('on');
                            
                    case 'annealing'
                        
                        nelderopts = optimset('Display','none',...
                                        'PlotFcns',[],...
                                        'MaxFunEvals',n_fitIterations,...
                                        'MaxIter',n_fitIterations,...
                                        'TolFun',1e-60,...
                                        'TolX',1e-60);
                                
                        annealopts = struct(...
                                        'CoolSched',@(T) (.4*T),...
                                        'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
                                        'InitTemp',1,...
                                        'MaxConsRej',1000,...
                                        'MaxSuccess',20,...
                                        'MaxTries',300,...
                                        'StopTemp',1e-5,...
                                        'StopVal',-Inf,...
                                        'Verbosity',0);
                        
                        if i == 1 && elasticSetting
                            % Fit the elastic term separately for the first
                            % iteration. Future iterations have the "best
                            % fit" elastic term included from the prior
                            % optimization attempt
                            
                            % Clock the timer
                            preElasticTime = toc;
                            
                            progressString = sprintf('Simulated Annealing (%s)\nInvestigating Elastic Parameter\nParallel Search Running...',model);
                            hbar = parfor_progressbar(n_iterations,progressString);
                            warning('off');
                            parfor j = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fminsearch(@(x)objFunc(x,elasticSetting,fluidSetting),beta0,nelderopts);
                                hbar.iterate(1);
                            end
                            close(hbar);
                            warning('on');
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = horzcat(fitStruct.elasticFitTime,{postElasticTime-preElasticTime});
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        progressString = sprintf('Simulated Annealing (%s)\nInvestigating Material Parameters\nParallel Search Running...',model);
                        hbar = parfor_progressbar(n_iterations,progressString);
                        warning('off');
                        parfor j = 1:n_iterations
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = annealOpt(@(x)objFunc(x,elasticSetting,fluidSetting),beta0,annealopts,nelderopts);
                            hbar.iterate(1);
                        end
                        close(hbar);
                        warning('on');
                        
                    case 'nls'
                        
                        fminoptions = optimoptions('fmincon','Algorithm','sqp',...
                                        'MaxFunctionEvaluations', n_fitIterations,...
                                        'MaxIterations', n_fitIterations,...
                                        'FiniteDifferenceType','central',...
                                        'FunctionTolerance', 1e-60,...
                                        'OptimalityTolerance', 1e-60,...
                                        'Display', 'none');
                        
                        if i == 1 && elasticSetting
                            % Fit the elastic term separately for the first
                            % iteration. Future iterations have the "best
                            % fit" elastic term included from the prior
                            % optimization attempt
                            
                            % Clock the timer
                            preElasticTime = toc;
                            
                            progressString = sprintf('Nonlinear Least-Squares (%s)\nInvestigating Elastic Parameter\nParallel Search Running...',model);
                            hbar = parfor_progressbar(n_iterations,progressString);
                            warning('off');
                            parfor j = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fmincon(@(x)objFunc(x,elasticSetting,fluidSetting),beta0,[],[],[],[],lb(1),ub(1),[],fminoptions);
                                hbar.iterate(1);
                            end
                            close(hbar);
                            warning('on');
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = horzcat(fitStruct.elasticFitTime,{postElasticTime-preElasticTime});
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        progressString = sprintf('Nonlinear Least-Squares (%s)\nInvestigating Material Parameters\nParallel Search Running...',model);
                        hbar = parfor_progressbar(n_iterations,progressString);
                        warning('off');
                        parfor j = 1:n_iterations
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = fmincon(@(x)objFunc(x,elasticSetting,fluidSetting),beta0,[],[],[],[],lb,ub,[],fminoptions);
                            hbar.iterate(1);
                        end
                        close(hbar);
                        warning('on');
                        
                    otherwise
                        error('That solver is not supported.')
                end
                postFitting = toc;
                
                % Find the best-fit parameters from our population
                [~,idx] = min(residual_dist,[],'omitnan');
                if size(idx,2)>1 idx = (idx(1)); end
                
                % Store the best-fit parameters to be fed forward, and also
                % save the "less optimal" sets for statistical treatment
                % later
                fprintf('Optimal Parameters, Iteration %d:\n',i)
                disp(beta_dist(:,idx));
                fitStruct.bestParams = horzcat(fitStruct.bestParams,{beta_dist(:,idx)});
                fitStruct.paramPopulation = horzcat(fitStruct.paramPopulation,{beta_dist});
                fitStruct.paramPopulationResiduals = horzcat(fitStruct.paramPopulationResiduals,{residual_dist});
                
                % Store the timing for this model configuration fit
                fitStruct.fitTime = horzcat(fitStruct.fitTime,{postFitting-preFitting});
                
                if any(strcmp(model,'plr'))
                    % The PLR model does not utilize more than a single
                    % term, so "iterative term introduction" has no
                    % meaning in this case. As such, the remaining
                    % iterations are skipped and the results are returned
                    
                    % To stop the iterative term introduction for one of
                    % the models, add an additional strcmp(model,'{name}')
                    % after the one in this if statement declaration,
                    % separated by a comma. The any() will trigger if the
                    % type is PLR or whatever new entry you have added.
                    % This is primarily in the case where you have a
                    % custom model that does not use iterative term
                    % introduction.
                    
                    break;
                end
                
            end % End Iterative Term Introduction Loop
            
            delete(poolobj);
            
        end % End fitData()
        
        function fitStruct = fitMap(obj, varargin)
            % FITMAP Fit a Viscoelastic Model to Force Map Data
            %   This function takes in a variety of settings in addition to
            %   the data already provided to the class and performs an
            %   optimization procedure based on those settings.
            %
            %   This function will result in an optimized parameter set
            %   for each pixel of a force map. To do this, the pixel
            %   treatment is parallelized such that each worker handles the
            %   fitting procedure for a single force curve. As opposed to
            %   having a single set of parameters for each model
            %   configuration, as occurs for fitData(), each model
            %   configuration will output a grid of parameter sets.
            %
            %   Defaults are below, and their order of definition
            %   must be adhered to (unless the user manually reprograms
            %   this class).
            %   
            %   In particular, the number of initializations chosen
            %   (n_iterations) may need to be increased, as with the number
            %   of solver iterations (n_fitIterations) and maximum number
            %   of elements in the viscoelastic series (n_elements). This
            %   function will iteratively introduce viscoelastic elements,
            %   feeding forward the results from the previous iteration,
            %   until the optimization has been performed for a model
            %   configuration with a number of viscoelastic elements 
            %   equal to n_elements.
            
            % Initialize Output Structure
            fitStruct = struct;
            
            % Default Settings
            solver = 'nelder-mead';     % Fit using Nelder-Mead Simplex
            model = 'maxwell';          % Use Generalized Maxwell Model
            n_elements = 3;             % Fit iteratively for up to 3 elements
            elasticSetting = 1;         % Include Elastic Term
            fluidSetting = 0;           % No Steady-State Fluidity
            n_iterations = 100;         % Use 100 random initializations
            n_fitIterations = 1e4;      % No. of iterations for solver
            if ~isempty(varargin)
                % Only one varargin is accepted, and it is a structure
                % containing all of the settings information we require
                fitOpts = varargin{1};
                
                % Get the settings provided
                solver = lower(fitOpts.solver);
                model = lower(fitOpts.model);
                n_elements = fitOpts.n_elements;
                elasticSetting = fitOpts.elasticSetting;
                fluidSetting = fitOpts.fluidSetting;
                n_iterations = fitOpts.n_iterations;
                n_fitIterations = fitOpts.n_fitIterations;
                if strcmpi(solver,'custom')
                    if isfield(fitOpts,'customFunc')
                        customFunc = fitOpts.customFunc;
                    else
                        error('You selected "custom" for the solver type and did not supply the function name in the class settings.');
                    end
                end
            end
            
            % Store the fit settings for future reference
            fitStruct.solver = solver;
            fitStruct.model = model;
            fitStruct.n_elements = n_elements;
            fitStruct.elasticSetting = elasticSetting;
            fitStruct.fluidSetting = fluidSetting;
            fitStruct.n_iterations = n_iterations;
            
            % Get the correct objective function for optimization
            switch model
                case 'maxwell'
                    objFuncMap = @obj.SSE_Maxwell_Map;
                    modelFile = 'LR_Maxwell.m';
                case 'voigt'
                    objFuncMap = @obj.SSE_Voigt_Map;
                    modelFile = 'LR_Voigt.m';
                case 'plr'
                    objFuncMap = @obj.SSE_PLR_Map;
                    modelFile = 'LR_PLR.m';
                    % The storage roster for plr is different, and requires
                    % the second index to maintain consistency. Thus, we
                    % have to manually force the second term to be fit by
                    % including the fluidity. This second position actually
                    % corresponds to the exponent, alpha.
                    fluidSetting = 1;
                    fitStruct.fluidSetting = fluidSetting;
                case 'custom'
                    objFuncMap = @obj.customFunc_Map;
                    modelFile = 'customFunc.m';
                otherwise
                    error('Your chosen solver-model combination is not implemented yet.');
            end
            
            % Create placeholders for the data we will obtain from our
            % optimization attempts. These need to be pre-allocated to
            % allow proper parallelization of the pixels.
            fitStruct.bestParams = cell(1,n_elements);
            fitStruct.bestParams(:) = cell(size(obj.forces_cell));
            
            fitStruct.paramPopulation = cell(1,n_elements);
            fitStruct.paramPopulation(:) = cell(size(obj.forces_cell));
            
            fitStruct.paramPopulationResiduals = cell(1,n_elements);
            fitStruct.paramPopulationResiduals(:) = cell(size(obj.forces_cell));
            
            fitStruct.elasticFitTime = cell(1,n_elements);
            fitStruct.elasticFitTime(:) = cell(size(obj.forces_cell));
            
            fitStruct.fitTime = cell(1,n_elements);
            fitStruct.fitTime(:) = cell(size(obj.forces_cell));
            
            % Open Parallel Pool of MATLAB Workers
            if isempty(gcp('nocreate'))
                % Make a fresh pool
                poolobj = parpool('IdleTimeout', Inf);
            else
                % Remove the old pool
                poolobj = gcp('nocreate');
                delete(poolobj);
                
                % Make a fresh pool
                poolobj = parpool('IdleTimeout', Inf);
            end
            
            % Send the class to the workers
            addAttachedFiles(poolobj, {'ViscoFit.m',modelFile})

            % Start the timer
            tic;
            
            % For Matlab's Parfor, we have to explicitly define the loop
            % bounds ahead of time:
            n_pixels = numel(obj.forces_cell);
            
            % We initialize the map variables that we will be updating. For
            % the first iteration, these are all initially empty.
            % Afterward, they will contain the previous model
            % configuration's results. These are updated in turn and stored
            % before being overwritten.
            bestParamsMap = cell(size(obj.forces_cell));
            paramPopulationMap = cell(size(obj.forces_cell));
            paramPopulationResidualsMap = cell(size(obj.forces_cell));
            elasticFitTimeMap = cell(size(obj.forces_cell));
            fitTimeMap = cell(size(obj.forces_cell));
            
            % Begin the iterative term introduction loop
            for i = 1:n_elements

                progressString = sprintf('Viscoelastic Force Map Analysis\n%d Viscoelastic Elements (%s config.)\nAnalyzing Pixels...',i,model);
                hbar = parfor_progressbar(n_pixels,progressString);
                warning('off');
                parfor j = 1:n_pixels

                    % Look to see if there are old results available to provide
                    % intelligent guesses for our initial conditions. This will
                    % only occur for iterations beyond the first
                    if i > 1
                        % Our current loop has an array of parameters
                        % two-larger than the previous (for the generalized
                        % spring-dashpot rheology models)
                        beta_in = NaN(length(bestParamsMap{j})+2,1);

                        % Overwrite the NaN values with the previous optimal
                        % parameters.
                        beta_in(1:length(bestParamsMap{j})) = bestParamsMap{j};
                    else
                        % For the first iteration, there will be four
                        % parameters in total: the elastic element (1), the
                        % steady-state fluidity (2), and two parameters for the
                        % first viscoelastic element (3,4). Note, based on the
                        % elasticSetting and fluidSetting provided to this
                        % function, the values in (1) and (2) may never be
                        % updated
                        if strcmpi(model,'plr')
                            beta_in = NaN(2,1);
                        else
                            beta_in = NaN(4,1);
                        end
                    end

                    % Make our storage variables
                    beta_dist = zeros(length(beta_in),n_iterations);
                    beta_dist_elastic = zeros(1,n_iterations);
                    beta0_dist = zeros(size(beta_dist));
                    residual_dist = NaN(1,n_iterations);
                    residual_dist_elastic = NaN(size(residual_dist));

                    % Make the upper and lower bounds for this model
                    [tauInds,modulusInds] = getParamIndices(beta_in);
                    ub = zeros(size(beta_in))+eps;
                    lb = zeros(size(beta_in));
                    switch model
                        case 'maxwell'
                            % Moduli are limited to "reasonable" bounds for
                            % viscoelastic materials
                            ub(modulusInds) = 1e12;
                            lb(modulusInds) = 1e-2;

                            tauCenters = obj.minTimescale.*(10.^( (1:length(ub(3:2:end)))-1 ));
                            ub(tauInds) = tauCenters*10;
                            lb(tauInds) = tauCenters/10;

                            if ~elasticSetting
                                ub(1) = eps;
                                lb(1) = 0;
                            end

                            if fluidSetting
                                ub(2) = max(tauCenters)*1e2;
                                lb(2) = 10^( floor(min(log10(obj.dts)))+1 );
                            end

                            % Restrict the range of random guesses, if desired.
                            % Otherwise, they should be set equal to ub & lb
                            ub_rand(modulusInds) = 6;
                            ub_rand(tauInds) = log10(ub(tauInds));
                            lb_rand(modulusInds) = 1;
                            lb_rand(tauInds) = log10(lb(tauInds));

                        case 'voigt'
                            % Compliances are limited to "reasonable" bounds
                            % for viscoelastic materials
                            ub(modulusInds) = 1e2;
                            lb(modulusInds) = 1e-12;

                            tauCenters = obj.minTimescale.*(10.^( (1:length(ub(3:2:end)))-1 ));
                            ub(tauInds) = tauCenters*10;
                            lb(tauInds) = tauCenters/10;

                            if ~elasticSetting
                                ub(1) = eps;
                                lb(1) = 0;
                            end

                            if fluidSetting
                                ub(2) = 1;
                                lb(2) = 0;
                            end

                            % Restrict the range of random guesses, if desired.
                            % Otherwise, they should be set equal to ub & lb
                            ub_rand(modulusInds) = -1;
                            ub_rand(tauInds) = log10(ub(tauInds));
                            lb_rand(modulusInds) = -6;
                            lb_rand(tauInds) = log10(lb(tauInds));

                        case 'plr'
                            % Power Law Rheology Roster:
                            % [E_0 alpha]
                            ub = [1e12 1];
                            lb = [1e-2 0];

                            % Restrict the range of random guesses, if desired.
                            % Otherwise, they should be set equal to ub & lb
                            ub_rand = log10(ub);
                            lb_rand = [log10(lb(1)) -3];

                        case 'custom'

                            % Define the upper and lower bounds for your custom
                            % function here. The output from this region should
                            % be two arrays, ub and lb, which are the same
                            % length as the number of parameters in the model.
                            % ub = [...];
                            % lb = [...];

                    end

                    preFitting = toc;
                    switch solver
                        case 'nelder-mead'

                            options = optimset('Display','none',...
                                        'PlotFcns',[],...
                                        'MaxFunEvals',n_fitIterations,...
                                        'MaxIter',n_fitIterations,...
                                        'TolFun',1e-60,...
                                        'TolX',1e-60);

                            if i == 1 && elasticSetting
                                % Fit the elastic term separately for the first
                                % iteration. Future iterations have the "best
                                % fit" elastic term included from the prior
                                % optimization attempt

                                % Clock the timer
                                preElasticTime = toc;

                                for k = 1:n_iterations
                                    % Get the grid search starting position
                                    beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{k});
                                    [beta_dist_elastic(k),residual_dist_elastic(k)] = fminsearch(@(x)objFuncMap(x,j,elasticSetting,fluidSetting),beta0,options);
                                    hbar.iterate(1);
                                end

                                % Clock the timer and save the fitting time
                                postElasticTime = toc;
                                elasticFitTimeMap{j} = postElasticTime-preElasticTime;

                                % Find the best elastic parameter
                                [~,idx] = min(residual_dist_elastic,[],'omitnan');
                                beta_in(1) = beta_dist_elastic(:,idx);
                            end

                            % See which parameters are new this time, so that
                            % information can be fed to our
                            % random-guess-generation function
                            newInds = isnan(beta_in);
                            
                            for k = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                                beta0_dist(:,k) = beta0;
                                [beta_dist(:,k),residual_dist(k)] = fminsearch(@(x)objFuncMap(x,j,elasticSetting,fluidSetting),beta0,options);
                                hbar.iterate(1);
                            end

                        case 'annealing'

                            nelderopts = optimset('Display','none',...
                                            'PlotFcns',[],...
                                            'MaxFunEvals',n_fitIterations,...
                                            'MaxIter',n_fitIterations,...
                                            'TolFun',1e-60,...
                                            'TolX',1e-60);

                            annealopts = struct(...
                                            'CoolSched',@(T) (.4*T),...
                                            'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
                                            'InitTemp',1,...
                                            'MaxConsRej',1000,...
                                            'MaxSuccess',20,...
                                            'MaxTries',300,...
                                            'StopTemp',1e-5,...
                                            'StopVal',-Inf,...
                                            'Verbosity',0);

                            if i == 1 && elasticSetting
                                % Fit the elastic term separately for the first
                                % iteration. Future iterations have the "best
                                % fit" elastic term included from the prior
                                % optimization attempt

                                % Clock the timer
                                preElasticTime = toc;

                                for k = 1:n_iterations
                                    % Get the grid search starting position
                                    beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{k});
                                    [beta_dist_elastic(k),residual_dist_elastic(k)] = fminsearch(@(x)objFuncMap(x,j,elasticSetting,fluidSetting),beta0,nelderopts);
                                    hbar.iterate(1);
                                end

                                % Clock the timer and save the fitting time
                                postElasticTime = toc;
                                elasticFitTimeMap{j} = postElasticTime-preElasticTime;

                                % Find the best elastic parameter
                                [~,idx] = min(residual_dist_elastic,[],'omitnan');
                                beta_in(1) = beta_dist_elastic(:,idx);
                            end

                            % See which parameters are new this time, so that
                            % information can be fed to our
                            % random-guess-generation function
                            newInds = isnan(beta_in);
                            
                            for k = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                                beta0_dist(:,k) = beta0;
                                [beta_dist(:,k),residual_dist(k)] = annealOpt(@(x)objFuncMap(x,j,elasticSetting,fluidSetting),beta0,annealopts,nelderopts);
                                hbar.iterate(1);
                            end

                        case 'nls'

                            fminoptions = optimoptions('fmincon','Algorithm','sqp',...
                                            'MaxFunctionEvaluations', n_fitIterations,...
                                            'MaxIterations', n_fitIterations,...
                                            'FiniteDifferenceType','central',...
                                            'FunctionTolerance', 1e-60,...
                                            'OptimalityTolerance', 1e-60,...
                                            'Display', 'none');

                            if i == 1 && elasticSetting
                                % Fit the elastic term separately for the first
                                % iteration. Future iterations have the "best
                                % fit" elastic term included from the prior
                                % optimization attempt

                                % Clock the timer
                                preElasticTime = toc;

                                for k = 1:n_iterations
                                    % Get the grid search starting position
                                    beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{k});
                                    [beta_dist_elastic(k),residual_dist_elastic(k)] = fmincon(@(x)objFuncMap(x,elasticSetting,fluidSetting),beta0,[],[],[],[],lb(1),ub(1),[],fminoptions);
                                    hbar.iterate(1);
                                end

                                % Clock the timer and save the fitting time
                                postElasticTime = toc;
                                elasticFitTimeMap{j} = postElasticTime-preElasticTime;

                                % Find the best elastic parameter
                                [~,idx] = min(residual_dist_elastic,[],'omitnan');
                                beta_in(1) = beta_dist_elastic(:,idx);
                            end

                            % See which parameters are new this time, so that
                            % information can be fed to our
                            % random-guess-generation function
                            newInds = isnan(beta_in);

                            for k = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                                beta0_dist(:,k) = beta0;
                                [beta_dist(:,k),residual_dist(k)] = fmincon(@(x)objFuncMap(x,elasticSetting,fluidSetting),beta0,[],[],[],[],lb,ub,[],fminoptions);
                                hbar.iterate(1);
                            end

                        otherwise
                            error('That solver is not supported.')
                    end
                    postFitting = toc;

                    % Find the best-fit parameters from our population
                    [~,idx] = min(residual_dist,[],'omitnan');
                    if size(idx,2)>1 idx = (idx(1)); end

                    % Store the best-fit parameters to be fed forward, and also
                    % save the "less optimal" sets for statistical treatment
                    % later
                    bestParamsMap{j} = beta_dist(:,idx);
                    paramPopulationMap{j} = beta_dist;
                    paramPopulationResidualsMap{j} = residual_dist;

                    % Store the timing for this model configuration fit
                    fitTimeMap{j} = postFitting-preFitting;

                end % End Pixel Loop
                close(hbar);
                warning('on');
                
                if any(strcmp(model,'plr'))
                    % The PLR model does not utilize more than a single
                    % term, so "iterative term introduction" has no
                    % meaning in this case. As such, the remaining
                    % iterations are skipped and the results are returned

                    % To stop the iterative term introduction for one of
                    % the models, add an additional strcmp(model,'{name}')
                    % after the one in this if statement declaration,
                    % separated by a comma. The any() will trigger if the
                    % type is PLR or whatever new entry you have added.
                    % This is primarily in the case where you have a
                    % custom model that does not use iterative term
                    % introduction.

                    break;
                end
                
                % Store updated values in output structure
                fitStruct.bestParams{i} = bestParamsMap;
                fitStruct.paramPopulation{i} = paramPopulationMap;
                fitStruct.paramPopulationResiduals{i} = paramPopulationResidualsMap;
                fitStruct.elasticFitTime{i} = elasticFitTimeMap;
                fitStruct.fitTime{i} = fitTimeMap;
            
            end % End Iterative Term Introduction Loop
                
            delete(poolobj);
            
        end % End fitMap()
        
    end % End Methods
    
end

