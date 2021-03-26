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
        times double {mustBePositive}
        dts double {mustBePositive}
        indentations double
        tipSize double {mustBePositive}
        tipGeom string
        minTimescale double {mustBePositive}
        nu double {mustBePositive}
    end
    
    methods
        function obj = ViscoFit(forces,times,indentations,tipSize,minTimescale,varargin)
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
            if ~isempty(varargin)
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
        end % End Maxwell SSE Function
        
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
        end % End Voigt SSE Function
        
        function sse = SSE_PLR(obj,params)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            
            % Calculate test forces
            test_forces = LR_PLR(params,obj.times,obj.dts,obj.indentations,obj.tipSize,obj.nu,obj.tipGeom);
            
            % calculate global residual
            sse_global = sum((obj.forces-test_forces).^2);
            sse = sum(sse_global);
            
            % Power Law Rheology Roster:
            % [E_0 alpha]
            ub = [1e12 1];
            lb = [1e-2 0];
            
            if any(ub-params < 0) || any(params-lb < 0)
                sse = Inf;
            end
        end % End PLR SSE Function
        
        function fitStruct = fit(obj, varargin)
            % FIT Fit a Viscoelastic Model to the Class Data 
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
            n_iterations = 100;         % Use 100 random initializations
            n_fitIterations = 1e4;      % No. of iterations for solver
            if ~isempty(varargin)
                for i=1:length(varargin)
                    switch(i)
                        case 1
                            solver = char(varargin{i});
                        case 2
                            model = char(varargin{i});
                        case 3
                            n_elements = cell2mat(varargin{i});
                        case 4
                            elasticSetting = cell2mat(varargin{i});
                        case 5
                            fluidSetting = cell2mat(varargin{i});
                        case 6
                            n_iterations = cell2mat(varargin{i});
                        case 7
                            n_fitIterations = cell2mat(varargin{i});
                        case 8
                            customFunc = varargin{i};
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
                    objFunc = @SSE_Maxwell;
                case 'voigt'
                    objFunc = @SSE_Voigt;
                case 'PLR'
                    objFunc = @SSE_PLR;
                case 'custom'
                    objFunc = @customFunc;
                otherwise
                    error('Your chosen solver-model combination is not implemented yet.');
            end
            
            % Create placeholders for the data we will obtain from our
            % optimization attempts
            fitStruct.bestParams = {};
            fitStruct.paramPopulation = {};
            fitStruct.paramPopulationResiduals = {};
            fitStruct.elasticFitTime = [];
            fitStruct.fitTime = {};
            
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
                    beta_in = NaN(4,1);
                end
                
                % Make our storage variables
                beta_dist = zeros(length(beta_in),n_iterations);
                beta_dist_elastic = zeros(1,size(beta_dist));
                beta0_dist = zeros(size(beta_dist));
                residual_dist = NaN(size(beta_dist));
                residual_dist_elastic = NaN(size(beta_dist_elastic));
                
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
                        ub_rand(modulusInds) = 10^(6);
                        ub_rand(tauInds) = ub(tauInds);
                        lb_rand(modulusInds) = 10^(1);
                        lb_rand(tauInds) = lb(tauInds);
                        
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
                        ub_rand(modulusInds) = 10^(-1);
                        ub_rand(tauInds) = ub(tauInds);
                        lb_rand(modulusInds) = 10^(-6);
                        lb_rand(tauInds) = lb(tauInds);
                        
                    case 'PLR'
                        % Power Law Rheology Roster:
                        % [E_0 alpha]
                        ub = [1e12 1];
                        lb = [1e-2 0];
                        
                        % Restrict the range of random guesses, if desired.
                        % Otherwise, they should be set equal to ub & lb
                        ub_rand = ub;
                        lb_rand = lb;
                        
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
                            hbar = parfor_progressbar(loopLim,progressString);
                            warning('off');
                            parfor j = 1:n_samples
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub(1),lb(1),n_samples),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fminsearch(@(x)objFunc(x),beta0,options);
                                hbar.iterate(1);
                            end
                            close(hbar);
                            warning('on');
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = postElasticTime-preElasticTime;
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        progressString = sprintf('Nelder-Mead (%s)\nInvestigating Material Parameters\nParallel Search Running...',model);
                        hbar = parfor_progressbar(loopLim,progressString);
                        warning('off');
                        parfor j = 1:n_samples
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = fminsearch(@(x)objFunc(x),beta0,options);
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
                                        'CoolSched',@(T) (.8*T),...
                                        'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
                                        'InitTemp',1,...
                                        'MaxConsRej',1000,...
                                        'MaxSuccess',20,...
                                        'MaxTries',n_fitIterations,...
                                        'StopTemp',scaling,...
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
                            hbar = parfor_progressbar(loopLim,progressString);
                            warning('off');
                            parfor j = 1:n_samples
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub(1),lb(1),n_samples),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fminsearch(@(x)objFunc(x),beta0,nelderopts);
                                hbar.iterate(1);
                            end
                            close(hbar);
                            warning('on');
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = postElasticTime-preElasticTime;
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        progressString = sprintf('Simulated Annealing (%s)\nInvestigating Material Parameters\nParallel Search Running...',model);
                        hbar = parfor_progressbar(loopLim,progressString);
                        warning('off');
                        parfor j = 1:n_samples
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = annealOpt(@(x)objFunc(x),beta0,annealopts,nelderopts);
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
                            hbar = parfor_progressbar(loopLim,progressString);
                            warning('off');
                            parfor j = 1:n_samples
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub(1),lb(1),n_samples),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fmincon(objFunc,beta0,[],[],[],[],lb(1),ub(1),[],fminoptions);
                                hbar.iterate(1);
                            end
                            close(hbar);
                            warning('on');
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = postElasticTime-preElasticTime;
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        progressString = sprintf('Nonlinear Least-Squares (%s)\nInvestigating Material Parameters\nParallel Search Running...',model);
                        hbar = parfor_progressbar(loopLim,progressString);
                        warning('off');
                        parfor j = 1:n_samples
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = fmincon(objFunc,beta0,[],[],[],[],lb(1),ub(1),[],fminoptions);
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
                
                % Store the best-fit parameters to be fed forward, and also
                % save the "less optimal" sets for statistical treatment
                % later
                fitStruct.bestParams = horzcat(fitStruct.bestParams,{beta_dist(:,idx)});
                fitStruct.paramPopulation = horzcat(fitStruct.paramPopulation,{beta_dist});
                fitStruct.paramPopulationResiduals = horzcat(fitStruct.paramPopulationResiduals,{residual_dist});
                
                % Store the timing for this model configuration fit
                fitStruct.fitTime = horzcat(fitStruct.fitTime,{postFitting-preFitting});
                
                if any(strcmp(model,'PLR'))
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
            
        end % End Fit()
        
    end % End Methods
    
end

