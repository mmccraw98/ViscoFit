clear all
close all
clc

addpath(genpath([pwd '\lib']));

% Pick the AFM Data Directory and Choose Data Extraction Settings
originalPath = uigetdir(pwd,...
        'Select the Folder Containing Your AFM Files');

% Check to see if there are subdirectories
dirContents = dir(originalPath);
subFolders = dirContents([dirContents.isdir]);
subFolders(contains({subFolders.name}, {'.','..','Plots'})) = [];

% If the user provides a main directory with many subdirectories containing
% data, we should loop through all directories and analyze each in turn.
if ~isempty(subFolders)
    Folders = cell(1,length(subFolders));
    Folders = cellfun(@(root,sub)[root '\' sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
else
    Folders = {originalPath};
end

for i_dir = 1:length(Folders)
    
    % Start with a clean slate
    clearvars -except i_dir Folders originalPath
    close all
    clc
    fprintf('Analyzing Directory #%d of %d\n',i_dir,length(Folders));
    
    % Use the current 
    path = Folders{i_dir};
    
    % Optional: Prompt user for minimum central fitting timescale and the
    % indenter geometry
    % minTimescale = input('Please enter the minimum timescale to use for fitting (e.g. 1e-4): ');
    % tipOpts = {"spherical","conical"};
    % [indx,~] = listdlg('PromptString','Indenter Geometry',...
    %     'SelectionMode','single',...
    %     'ListString',tipOpts);
    % tipGeom = tipOpts{indx};
    
    % Alternatively, set those values manually
    minTimescale = 1e-4;                    % This is the time value on which the first viscoelastic element will be centered
    tipGeom = "spherical";                  % The experiment tip geometry for the files that are being loaded
    
    % Settings for how to use the loaded data during analysis
    useSmoothData = 0;                      % The user can choose to use filtered data, or the original raw data for fitting
    useAveragedData = 1;                    % Choose whether to average force curves with the same approach velocity 
    
    % Settings for loading the data
    loadDataSettings  = struct();
    
    % Required Settings:
    loadDataSettings.includeRetract = 1;             % Don't include data from the retract curve
    loadDataSettings.filterType = 'none';            % Choose the filter used to smooth data
    loadDataSettings.findRep = 'forward';            % Search direction for the repulsive region
    loadDataSettings.removeNegatives = 1;            % Remove negative values in the data stream
    
    % Conditional Settings (depending on filter):
    loadDataSettings.N = 2;                          % Order of Butterworth filter (if used)
    loadDataSettings.cutoff_Hz = 5000;               % Cutoff frequency of Butterworth (if used)
    
    % Load the AFM Data
    dataStruct = LoadAFMData(path,loadDataSettings);

    % Choose the right data from our structure
    % To begin, we need to either select the averaged data rows (stored at the
    % bottom of the structure), and know the number of files used to make them
    % (the rows leading up to the averages). There will be one averaged
    % datasets per approach velocity.
    numFiles = 0;
    avgCount = 0;
    for i = 1:size(dataStruct,2)
        if isempty(dataStruct(i).t_average)
            numFiles = numFiles + 1;
        else
            avgCount = avgCount + 1;
        end
    end

    % Get our offsets, knowing how many files and averages there are
    if useAveragedData
        indShift = numFiles;
        loopMax = avgCount;
    else
        indShift = 0;
        loopMax = numFiles;
    end

    % Create our cell arrays containing the data we care to send to the
    % ViscoFit class
    forces = cell(1,loopMax);
    times = cell(size(forces));
    indentations = cell(size(forces));
    tipSize = cell(size(forces));
    nu = cell(size(forces));

    for i = 1:loopMax
        if useSmoothData
            forces{i} = dataStruct(indShift+i).F_r_smooth;
            times{i} = dataStruct(indShift+i).t_r_smooth;
            indentations{i} = dataStruct(indShift+i).h_r_smooth;
        else
            forces{i} = dataStruct(indShift+i).F_r;
            times{i} = dataStruct(indShift+i).t_r;
            indentations{i} = dataStruct(indShift+i).h_r;
        end
        tipSize{i} = dataStruct(indShift+i).r_tip;
        nu{i} = dataStruct(indShift+i).nu_sample;
    end

    % Test the Fitting Functions using the Nelder-Mead Solver
    % Make a structure for our settings
    classSettings = struct;

    % Test the Maxwell
    classSettings.minTimescale = 1e-4;              % Timescale for the first viscoelastic element
    classSettings.nu = nu;                          % Sample Poisson Ratio for all curves
    classSettings.tipGeom = tipGeom;                % Tip geometry for these experiments
    classSettings.fitLog = false;                   % Log-scale-resample the data before fitting (faster)

    % Create the class object
    visco = ViscoFit(forces,times,indentations,tipSize,classSettings);

    % Make a structure for our settings
    fitSettings = struct;

    % Test the Maxwell
    fitSettings.solver = 'nelder-mead';             % Fit using Nelder-Mead Simplex
    fitSettings.model = 'maxwell';                  % Use Generalized Maxwell Model
    fitSettings.n_elements = 4;                     % Fit iteratively for up to 4 elements
    fitSettings.elasticSetting = 1;                 % Include Elastic Term
    fitSettings.fluidSetting = 0;                   % No Steady-State Fluidity
    fitSettings.n_iterations = 200;                 % Use 200 random initializations
    fitSettings.n_fitIterations = 5e3;              % No. of iterations for solver
    maxwellFit_NM = visco.fitData(fitSettings);

    % Test the Voigt
    fitSettings.model = 'voigt';
    voigtFit_NM = visco.fitData(fitSettings);

    % Test the PLR
    fitSettings.model = 'PLR';
    PLRFit_NM = visco.fitData(fitSettings);

    % For a custum function, use the below code to pass the appropriate
    % function call to the fitting class. Note that an error will be thrown if
    % you select "custom" for the model type and then don't pass a function
    % name in the settings.
    % fitSettings.model = 'custom';
    % fitSettings.customFunc = @customFuncName;
    % customFit = visco.fitData(fitSettings);

    % Check out the results
    dts = cellfun(@(t) round(mode(gradient(t)),1,'significant').*ones(size(t)),times,'UniformOutput',false);
    switch tipGeom
        case "spherical"
            beta = 1.5;
        case "conical"
            beta = 2;
    end
    if exist('resultsFigNelder','var')
        try
            close(resultsFigNelder);
        catch
        end
        clearvars resultsFigNelder
    end
    resultsFigNelder = figure('Position',[50 100 300*fitSettings.n_elements 600]);
    for i = 1:fitSettings.n_elements
        subplot(1,fitSettings.n_elements,i)
        title(sprintf('%d Terms',i))
        hold on
        for j = 1:size(forces,2)
            scatter(times{j},forces{j},50,'rx')
            scatter(times{j},indentations{j}.^(beta),50,'bo')
            plot(times{j},LR_Maxwell(maxwellFit_NM.bestParams{i},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'r-','linewidth',3)
            plot(times{j},LR_Voigt(voigtFit_NM.bestParams{i},times{j},dts{j},forces{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'b-','linewidth',3)
            plot(times{j},LR_PLR(PLRFit_NM.bestParams{1},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'g-','linewidth',3)
        end
        grid on
        set(gca,'xscale','log','yscale','log')
        hold off
    end
    saveas(resultsFigNelder,[path '\PlotResults-NelderMead.jpg']);
    saveas(resultsFigNelder,[path '\PlotResults-NelderMead.fig']);
    save([path '\FitResults-NelderMead.mat'],'maxwellFit_NM','voigtFit_NM','PLRFit_NM')

    % Test the Fitting Functions using Simulated Annealing with Nelder-Mead
    % Create the class object
    % visco = ViscoFit(forces,times,indentations,tipSize,minTimescale,nu,tipGeom);

    % Make a structure for our settings
    % fitSettings = struct;

    % Test the Maxwell
    fitSettings.solver = 'annealing';               % Fit using Simulated Annealing
    fitSettings.model = 'maxwell';                  % Use Generalized Maxwell Model
    fitSettings.n_elements = 4;                     % Fit iteratively for up to 4 elements
    fitSettings.elasticSetting = 1;                 % Include Elastic Term
    fitSettings.fluidSetting = 0;                   % No Steady-State Fluidity
    fitSettings.n_iterations = 5;                   % Use 5 random initializations
    fitSettings.n_fitIterations = 5e2;              % No. of iterations for solver
    maxwellFit_Anneal = visco.fitData(fitSettings);

    % Test the Voigt
    fitSettings.model = 'voigt';
    voigtFit_Anneal = visco.fitData(fitSettings);

    % Test the PLR
    fitSettings.model = 'PLR';
    PLRFit_Anneal = visco.fitData(fitSettings);

    % For a custum function, use the below code to pass the appropriate
    % function call to the fitting class. Note that an error will be thrown if
    % you select "custom" for the model type and then don't pass a function
    % name in the settings.
    % fitSettings.model = 'custom';
    % fitSettings.customFunc = @customFuncName;
    % customFit = visco.fitData(fitSettings);

    % Check out the results
    dts = cellfun(@(t) mode(gradient(t)).*ones(size(t)),times,'UniformOutput',false);
    switch tipGeom
        case "spherical"
            beta = 1.5;
        case "conical"
            beta = 2;
    end
    if exist('resultsFigAnnealing','var')
        try
            close(resultsFigAnnealing);
        catch
        end
        clearvars resultsFigAnnealing
    end
    resultsFigAnnealing = figure('Position',[50 200 300*fitSettings.n_elements 600]);
    for i = 1:fitSettings.n_elements
        subplot(1,fitSettings.n_elements,i)
        title(sprintf('%d Terms',i))
        hold on
        for j = 1:size(forces,2)
            scatter(times{j},forces{j},50,'rx')
            scatter(times{j},indentations{j}.^(beta),50,'bo')
            plot(times{j},LR_Maxwell(maxwellFit_Anneal.bestParams{i},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'r-','linewidth',3)
            plot(times{j},LR_Voigt(voigtFit_Anneal.bestParams{i},times{j},dts{j},forces{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'b-','linewidth',3)
            plot(times{j},LR_PLR(PLRFit_Anneal.bestParams{1},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'g-','linewidth',3)
        end
        grid on
        set(gca,'xscale','log','yscale','log')
        hold off
    end
    saveas(resultsFigAnnealing,[path '\PlotResults-Annealing.jpg']);
    saveas(resultsFigAnnealing,[path '\PlotResults-Annealing.fig']);
    save([path '\FitResults-Annealing.mat'],'maxwellFit_Anneal','voigtFit_Anneal','PLRFit_Anneal')

    % Test the Fitting Functions using Nonlinear Least Squares (lsqcurvefit)
    % Create the class object
    % visco = ViscoFit(forces,times,indentations,tipSize,minTimescale,nu,tipGeom);

    % Make a structure for our settings
    % fitSettings = struct;

    % Test the Maxwell
    fitSettings.solver = 'nls';                     % Fit using lsqcurvefit
    fitSettings.model = 'maxwell';                  % Use Generalized Maxwell Model
    fitSettings.n_elements = 4;                     % Fit iteratively for up to 4 elements
    fitSettings.elasticSetting = 1;                 % Include Elastic Term
    fitSettings.fluidSetting = 0;                   % No Steady-State Fluidity
    fitSettings.n_iterations = 200;                 % Use 200 random initializations
    fitSettings.n_fitIterations = 1e4;              % No. of iterations for solver
    maxwellFit_NLS = visco.fitData(fitSettings);

    % Test the Voigt
    fitSettings.model = 'voigt';
    voigtFit_NLS = visco.fitData(fitSettings);

    % Test the PLR
    fitSettings.model = 'PLR';
    PLRFit_NLS = visco.fitData(fitSettings);

    % For a custum function, use the below code to pass the appropriate
    % function call to the fitting class. Note that an error will be thrown if
    % you select "custom" for the model type and then don't pass a function
    % name in the settings.
    % fitSettings.model = 'custom';
    % fitSettings.customFunc = @customFuncName;
    % customFit = visco.fitData(fitSettings);

    % Check out the results
    dts = cellfun(@(t) mode(gradient(t)).*ones(size(t)),times,'UniformOutput',false);
    switch tipGeom
        case "spherical"
            beta = 1.5;
        case "conical"
            beta = 2;
    end
    if exist('resultsFigNLS','var')
        try
            close(resultsFigNLS);
        catch
        end
        clearvars resultsFigNLS
    end
    resultsFigNLS = figure('Position',[50 300 300*fitSettings.n_elements 600]);
    for i = 1:fitSettings.n_elements
        subplot(1,fitSettings.n_elements,i)
        title(sprintf('%d Terms',i))
        hold on
        for j = 1:size(forces,2)
            scatter(times{j},forces{j},50,'rx')
            scatter(times{j},indentations{j}.^(beta),50,'bo')
            plot(times{j},LR_Maxwell(maxwellFit_NLS.bestParams{i},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'r-','linewidth',3)
            plot(times{j},LR_Voigt(voigtFit_NLS.bestParams{i},times{j},dts{j},forces{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'b-','linewidth',3)
            plot(times{j},LR_PLR(PLRFit_NLS.bestParams{1},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'g-','linewidth',3)
        end
        grid on
        set(gca,'xscale','log','yscale','log')
        hold off
    end
    saveas(resultsFigNLS,[path '\PlotResults-NLS.jpg']);
    saveas(resultsFigNLS,[path '\PlotResults-NLS.fig']);
    save([path '\FitResults-NLS.mat'],'maxwellFit_NLS','voigtFit_NLS','PLRFit_NLS')
    
end

% Open the originally requested directory
winopen(originalPath);