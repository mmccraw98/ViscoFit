clear all
close all
clc

addpath(genpath([pwd '\lib']));

%% Pick the AFM Data Directory and Choose Data Extraction Settings
path = uigetdir(pwd,...
        'Select the Folder Containing Your AFM Files');
minTimescale = input('Please enter the minimum timescale to use for fitting (e.g. 1e-4): ');
tipOpts = {"spherical","conical"};
[indx,~] = listdlg('PromptString','Indenter Geometry',...
    'SelectionMode','single',...
    'ListString',tipOpts);
tipGeom = tipOpts{indx};
useSmoothData = 0;
useAveragedData = 1;

%% Load the AFM Data
dataStruct = LoadAFMData(path);

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

%% Test the Fitting Functions using the Nelder-Mead Solver
% Create the class object
visco = ViscoFit(forces,times,indentations,tipSize,minTimescale,nu,tipGeom);

% Make a structure for our settings
fitSettings = struct;

% Test the Maxwell
fitSettings.solver = 'nelder-mead';     % Fit using Nelder-Mead Simplex
fitSettings.model = 'maxwell';          % Use Generalized Maxwell Model
fitSettings.n_elements = 4;             % Fit iteratively for up to 4 elements
fitSettings.elasticSetting = 1;         % Include Elastic Term
fitSettings.fluidSetting = 0;           % No Steady-State Fluidity
fitSettings.n_iterations = 500;         % Use 500 random initializations
fitSettings.n_fitIterations = 1e4;      % No. of iterations for solver
maxwellFit = visco.fitData(fitSettings);

% Test the Voigt
fitSettings.model = 'voigt';
voigtFit = visco.fitData(fitSettings);

% Test the PLR
fitSettings.model = 'PLR';
PLRFit = visco.fitData(fitSettings);

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
if exist('resultsFig','var')
    try
        close(resultsFig);
    catch
    end
    clearvars resultsFig
end
resultsFig = figure('Position',[50 200 300*fitSettings.n_elements 600]);
for i = 1:fitSettings.n_elements
    subplot(1,fitSettings.n_elements,i)
    title(sprintf('%d Terms',i))
    hold on
    for j = 1:size(forces,2)
        scatter(times{j},forces{j},50,'rx')
        scatter(times{j},indentations{j}.^(beta),50,'bo')
        plot(times{j},LR_Maxwell(maxwellFit.bestParams{j},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'r-','linewidth',3)
        plot(times{j},LR_Voigt(voigtFit.bestParams{j},times{j},dts{j},forces{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'b-','linewidth',3)
        plot(times{j},LR_PLR(PLRFit.bestParams{1},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'g-','linewidth',3)
    end
    grid on
    set(gca,'xscale','log','yscale','log')
    hold off
end

%% Test the Fitting Functions using Simulated Annealing with Nelder-Mead
% Create the class object
% visco = ViscoFit(forces,times,indentations,tipSize,minTimescale,nu,tipGeom);

% Make a structure for our settings
% fitSettings = struct;

% Test the Maxwell
fitSettings.solver = 'annealing';       % Fit using Simulated Annealing
fitSettings.model = 'maxwell';          % Use Generalized Maxwell Model
fitSettings.n_elements = 4;             % Fit iteratively for up to 4 elements
fitSettings.elasticSetting = 1;         % Include Elastic Term
fitSettings.fluidSetting = 0;           % No Steady-State Fluidity
fitSettings.n_iterations = 20;          % Use 20 random initializations
fitSettings.n_fitIterations = 2.5e3;    % No. of iterations for solver
maxwellFit = visco.fitData(fitSettings);

% Test the Voigt
fitSettings.model = 'voigt';
voigtFit = visco.fitData(fitSettings);

% Test the PLR
fitSettings.model = 'PLR';
PLRFit = visco.fitData(fitSettings);

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
if exist('resultsFig','var')
    try
        close(resultsFig);
    catch
    end
    clearvars resultsFig
end
resultsFig = figure('Position',[50 200 300*fitSettings.n_elements 600]);
for i = 1:fitSettings.n_elements
    subplot(1,fitSettings.n_elements,i)
    title(sprintf('%d Terms',i))
    hold on
    for j = 1:size(forces,2)
        scatter(times{j},forces{j},50,'rx')
        scatter(times{j},indentations{j}.^(beta),50,'bo')
        plot(times{j},LR_Maxwell(maxwellFit.bestParams{j},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'r-','linewidth',3)
        plot(times{j},LR_Voigt(voigtFit.bestParams{j},times{j},dts{j},forces{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'b-','linewidth',3)
        plot(times{j},LR_PLR(PLRFit.bestParams{1},times{j},dts{j},indentations{j},tipSize{j},nu{j},tipGeom,fitSettings.elasticSetting,fitSettings.fluidSetting),'g-','linewidth',3)
    end
    grid on
    set(gca,'xscale','log','yscale','log')
    hold off
end