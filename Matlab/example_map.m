% Make a clean slate
clear all
close all
clc

% Add the library to our path
addpath(genpath([pwd filesep 'lib']));

% Pick the AFM Data Directory and Choose Data Extraction Settings
originalPath = uigetdir(pwd,...
        'Select the Folder Containing Your AFM Files');
savePrependIn = 'TESTMAP';

% Check to see if there are subdirectories
dirContents = dir(originalPath);
subFolders = dirContents([dirContents.isdir]);
subFolders(contains({subFolders.name}, {'.','..','Plots'})) = [];

% If the user provides a main directory with many subdirectories containing
% data, we should loop through all directories and analyze each in turn.
if ~isempty(subFolders)
    Folders = cell(1,length(subFolders));
    Folders = cellfun(@(root,sub)[root filesep sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
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
    tipGeom = "conical";                    % The experiment tip geometry for the files that are being loaded
    
    % Settings for how to use the loaded data during analysis
    useSmoothData = 0;                      % The user can choose to use filtered data, or the original raw data for fitting
    useAveragedData = 0;                    % Choose whether to average force curves with the same approach velocity 
    
    % Settings for loading the data
    loadDataSettings  = struct();
    
    % Required Settings:
    loadDataSettings.includeRetract = 0;             % Don't include data from the retract curve
    loadDataSettings.filterType = 'FIR';             % Choose the filter used to smooth data
    loadDataSettings.findRep = 'reverse';            % Search direction for the repulsive region
    loadDataSettings.removeNegatives = true;         % Remove negative values in the data stream
    loadDataSettings.createAverage = false;          % Create averaged rows at the END of the datastruct
    
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
        if ~isempty(dataStruct(i).z)
            numFiles = numFiles + 1;
        else
            avgCount = avgCount + 1;
        end
    end

    % Get our offsets, knowing how many files and averages there are
    if useAveragedData
        if avgCount == 0
            error('There are no averages in the dataStruct. You selected Averaged analysis. Please reload the data with "createAverage"==true.');
        end
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
    mapSize = cell(size(forces));

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
        mapSize{i} = dataStruct(indShift+i).mapSize;
    end
    
    % Memory Management
    clearvars dataStruct

    % Test the Fitting Functions using the Nelder-Mead Solver
    % Make a structure for our settings
    classSettings = struct;
    
    timescaleEst = 10^mode(ceil(log10(gradient([times{:}]))));

    % Test the Maxwell
    classSettings.minTimescale = 10*timescaleEst;   % Timescale for the first viscoelastic element
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
    fitSettings.n_elements = 1;                     % Fit iteratively for up to 3 elements
    fitSettings.elasticSetting = 1;                 % Include Elastic Term
    fitSettings.fluidSetting = 0;                   % No Steady-State Fluidity
    fitSettings.n_iterations = 10;                  % Use 10 random initializations
    fitSettings.n_fitIterations = 1e3;              % No. of iterations for solver
    fitSettings.errortype = 'mse';                  % Use Mean-Squared Error during fitting
    maxwellFit_NM = visco.fitMap(fitSettings);
    maxwellFit_NM.mapSize = mode(mapSize);
    save([path filesep savePrependIn '-FitResults-NelderMead.mat'],'maxwellFit_NM','-v7.3')
    clearvars maxwellFit_NM
    
    % Test the Voigt
    fitSettings.model = 'voigt';
    voigtFit_NM = visco.fitMap(fitSettings);
    voigtFit_NM.mapSize = mode(mapSize);
    save([path filesep savePrependIn '-FitResults-NelderMead.mat'],'voigtFit_NM','-append')
    clearvars voigtFit_NM

    % Test the PLR
    fitSettings.model = 'PLR';
    PLRFit_NM = visco.fitMap(fitSettings);
    PLRFit_NM.mapSize = mode(mapSize);
    save([path filesep savePrependIn '-FitResults-NelderMead.mat'],'PLRFit_NM','-append')
    clearvars PLRFit_NM
    
    % Make a pointer for our now MASSIVE results file
    outMat = matfile([path filesep savePrependIn '-FitResults-NelderMead.mat']);

    % For a custum function, use the below code to pass the appropriate
    % function call to the fitting class. Note that an error will be thrown if
    % you select "custom" for the model type and then don't pass a function
    % name in the settings.
    % fitSettings.model = 'custom';
    % fitSettings.customFunc = @customFuncName;
    % customFit = visco.fitMap(fitSettings);

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
    resultsFigNelder = figure('Position',[50 100 1200 600]);
    maxwellParams = outMat.maxwellFit_NM.bestParams;
    voigtParams = outMat.voigtFit_NM.bestParams;
    PLRParams = outMat.PLRFit_NM.bestParams;
    for i = 1:fitSettings.n_elements
        for j = 1:size(forces,2)
            scatter(times{j},forces{j},50,'rx')
            scatter(times{j},indentations{j}.^(beta),50,'bo')
            subplot(1,3,1)
            title('GM')
            hold on
            plot(times{j},LR_Maxwell(maxwellParams{i}{j},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'r-','linewidth',3)
            hold off
            subplot(1,3,2)
            title('GKV')
            hold on
            plot(times{j},LR_Voigt(voigtParams{i}{j},times{j},dts{j},forces{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'b-','linewidth',3)
            hold off
            subplot(1,3,3)
            title('PLR')
            hold on
            plot(times{j},LR_PLR(PLRParams{1}{j},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'g-','linewidth',3)
            hold off
        end
        grid on
        set(gca,'xscale','log','yscale','log')
        hold off
        saveas(resultsFigNelder,[path filesep savePrependIn sprintf('-PlotResults-NelderMead-%d_terms.jpg',i)]);
        saveas(resultsFigNelder,[path filesep savePrependIn sprintf('-PlotResults-NelderMead-%d_terms.fig',i)]);
        clf(resultsFigNelder)
    end
    clearvars maxwellParams voigtParams PLRParams outMat
    close(resultsFigNelder)

%     % Test the Fitting Functions using Nonlinear Least Squares (lsqcurvefit)
%     % Create the class object
%     % visco = ViscoFit(forces,times,indentations,tipSize,minTimescale,nu,tipGeom);
% 
%     % Make a structure for our settings
%     % fitSettings = struct;
% 
%     % Test the Maxwell
%     fitSettings.solver = 'nls';                     % Fit using lsqcurvefit
%     fitSettings.model = 'maxwell';                  % Use Generalized Maxwell Model
%     fitSettings.n_iterations = 200;                 % Use 200 random initializations
%     fitSettings.n_fitIterations = 1e4;              % No. of iterations for solver
%     maxwellFit_NLS = visco.fitMap(fitSettings);
% 
%     % Test the Voigt
%     fitSettings.model = 'voigt';
%     voigtFit_NLS = visco.fitMap(fitSettings);
% 
%     % Test the PLR
%     fitSettings.model = 'PLR';
%     PLRFit_NLS = visco.fitMap(fitSettings);
% 
%     % For a custum function, use the below code to pass the appropriate
%     % function call to the fitting class. Note that an error will be thrown if
%     % you select "custom" for the model type and then don't pass a function
%     % name in the settings.
%     % fitSettings.model = 'custom';
%     % fitSettings.customFunc = @customFuncName;
%     % customFit = visco.fitMap(fitSettings);
% 
%     % Check out the results
%     dts = cellfun(@(t) mode(gradient(t)).*ones(size(t)),times,'UniformOutput',false);
%     switch tipGeom
%         case "spherical"
%             beta = 1.5;
%         case "conical"
%             beta = 2;
%     end
%     if exist('resultsFigNLS','var')
%         try
%             close(resultsFigNLS);
%         catch
%         end
%         clearvars resultsFigNLS
%     end
%     resultsFigNLS = figure('Position',[50 300 300*fitSettings.n_elements 600]);
%     for i = 1:fitSettings.n_elements
%         subplot(1,fitSettings.n_elements,i)
%         title(sprintf('%d Terms',i))
%         hold on
%         for j = 1:size(forces,2)
%             scatter(times{j},forces{j},50,'rx')
%             scatter(times{j},indentations{j}.^(beta),50,'bo')
%             plot(times{j},LR_Maxwell(maxwellFit_NLS.bestParams{i},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'r-','linewidth',3)
%             plot(times{j},LR_Voigt(voigtFit_NLS.bestParams{i},times{j},dts{j},forces{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'b-','linewidth',3)
%             plot(times{j},LR_PLR(PLRFit_NLS.bestParams{1},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'g-','linewidth',3)
%         end
%         grid on
%         set(gca,'xscale','log','yscale','log')
%         hold off
%     end
%     saveas(resultsFigNLS,[path filesep 'PlotResults-NLS.jpg']);
%     saveas(resultsFigNLS,[path filesep 'PlotResults-NLS.fig']);
%     save([path filesep 'FitResults-NLS.mat'],'maxwellFit_NLS','voigtFit_NLS','PLRFit_NLS')
    
end

%% Open the originally requested directory
winopen(originalPath);
