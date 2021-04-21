clear all
close all
clc

addpath(genpath([pwd '\lib']));
addpath(genpath([pwd '\plotting']));

stillRunning = true;
i = 1;
while stillRunning
    
    if exist('pathList','var')
        startdir = pathList{i-1};
    else
        startdir = pwd;
    end
    
    % Pick the AFM Data Directory to put on the plot
    pathList{i} = uigetdir(startdir,...
            sprintf('Select Folder #%d to Compare',i));
    i = i + 1;
    
    % Prompt user
    answer = questdlg('Would you like to add another directory?', ...
        'Options', ...
        'No','Yes','Yes');

    % Handle response
    switch answer
        case 'No'
            clearvars -except pathList
            stillRunning = false;
        case 'Yes'
            close all
            clearvars -except stillRunning pathList i
    end
    
end

% Choose and output directory for the plots
outputpathname = uigetdir(pathList{end},...
        'Select the Output Directory for your Bar Plot');


%% Perform the Plotting Steps
% This section has been separated from the Data Loading section above so
% the user can easily run this section again after plotting to see the
% results for the SAME directories but different settings (e.g., solver
% choice).

% Choose which solver results to use
optList = {'Nelder-Mead','Annealing','NLS'}; 
[idx, tf] = listdlg('ListString', optList,...
    'SelectionMode', 'Single',...
    'PromptString', 'Choose which solver you used for the fit',...
    'Initialvalue', 1,...
    'Name', 'Select Your Solver');

if tf
    solverChoice = optList{idx};
else
    error('Please select which solver you used for the fit. This setting is required for this visualization script.');
end

% Clear old figures if they exist and make according to the below size
% settings
figX = 400;
figY = 50;
figWid = 1600;
figHeight = 800;
if ~exist('errorPlot','var')
    errorPlot = figure('Position',[figX figY figWid figHeight]);
else
    try
        figure(errorPlot)
        clf
    catch
        clearvars errorPlot
        errorPlot = figure('Position',[figX figY figWid figHeight]);
    end
end

% Decide how many subplots we will have based on the number of directories
n_subplots = length(pathList);

for i_paths = 1:length(pathList)
    
    % User-Defined Settings
    nSims = 36;                 % Number of simulations total
    nCols = 1;
    errortype = 'harmonic';     % harmonics or time
    errorfunc = 'mse';          % sse or mse
    linestyleList = {'-','--',':'};
    
    % Check to see if there are subdirectories
    dirContents = dir(pathList{i_paths});
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

    % Make a list of unique colors for the test conditions
    testConditionColors = num2cell(jet(nSims),2);

    % Pre-process some results to detect how many subplots we need
    nRows = ceil(n_subplots/nCols);
    
    % Define error functions
    sse_global = @(data,model) sum((data-model).^2,'all');
    mse_global = @(data,model,n) sum(((data-model).^2)./(movvar(data,3).^2),'all')./(length(data)-n);
        
    % Select the correct subplot (i.e. directory)
    subplot(nRows,nCols,i_paths)
    
    % Initialize Array for Bar Plot Data
    % There is one column for each model
    barData = NaN(0,3);
    numArmData = NaN(0,3);
    armLabelColorData = cell(0,3);
    testCondList = [];
    testCondData = {};
    
    % Begin looping through the directories or files
    for i_dir = 1:length(Folders)
        
        path = Folders{i_dir};
        Files = dir([path '\FitResults*.mat']);

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

        minTime = 0;
        maxTime = Inf;
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
            if min(times{i}) > minTime
                minTime = min(times{i});
            end
            if max(times{i}) < maxTime
                maxTime = max(times{i});
            end
        end

        % Generate a frequency array in log scale
        dt = mode(cell2mat(cellfun(@(x)mode(round(gradient(x),4,'significant')),times,'UniformOutput',false)));
        wmin = log10(2.0.*pi.*(1./maxTime));
        wmax = log10(2.0.*pi.*(1./dt));
        omega = logspace(wmin,wmax,round(10*(wmax-wmin)));

        clearvars dataStruct
        
        % Initialize results for this loop
        barDataLoop = NaN(1,3);
        numArmLoop = NaN(1,3);
        armLabelColorLoop = cell(1,3);
        
        % Get simulation parameters
        tempDir = strsplit(path,{'\','/'});
        tempStr = cell2mat(regexp(tempDir{end},'TestCondition[0-9]*','match'));
        simNum = str2double(cell2mat(regexp(tempStr,'\d*','match')));

        for j_dir = 1:length(Files)
            fitResults = load([Files(j_dir).folder '\' Files(j_dir).name],'-mat');
            dataLabels = fieldnames(fitResults);
            resultsContainer = cell(numel(dataLabels),1);
            for j = 1:numel(dataLabels)
                resultsContainer{j} = fitResults.(dataLabels{j});
            end

            simParams = getSimParams(simNum);
            simParams = horzcat(simParams(1),0,simParams(2:end));

            simSettings = struct;
            simSettings.model = 'maxwell';
            simSettings.bestParams = simParams;
            simSettings.elasticSetting = true;
            simSettings.fluidSetting = false;
            simSettings.n_elements = (numel(simParams)-2)/2;
                        
            for j = 1:numel(resultsContainer)
                tempResults = resultsContainer{j};
                
                if ~strcmpi(tempResults.solver,solverChoice)
                    continue;
                end

                tempSettings = struct;
                tempSettings.minTimescale = 1e-4;
                tempSettings.nu = nu;
                tempSettings.tipGeom = "spherical";
                tempSettings.fitLog = false;
                visco = ViscoFit(forces,times,indentations,tipSize,tempSettings);
                [simStorage,simLoss,simAng] = visco.harmonics_Maxwell(omega,simSettings);
                simAbsMod = sqrt(simStorage.^2+simLoss.^2);

                % Initialize Best Error Record (for number of model arms)
                loopErrorRecord = Inf(1,3);
                
                for k = 1:numel(tempResults.bestParams)

                    harmonicSettings = struct;
                    harmonicSettings.bestParams = tempResults.bestParams{k};
                    harmonicSettings.elasticSetting = tempResults.elasticSetting;
                    harmonicSettings.fluidSetting = tempResults.fluidSetting;
                    harmonicSettings.model = tempResults.model;

                    markerSize = 100+numel(tempResults.bestParams)*100*(1-k/numel(tempResults.bestParams));

                    switch lower(tempResults.model)
                        case 'maxwell'
                            recordColNum = 1;
                            [modelStorage,modelLoss,~] = visco.harmonics_Maxwell(omega,harmonicSettings);
                            switch errorfunc
                                case 'sse'
                                    modelErrorHarmonic = sse_global(simAbsMod,sqrt(modelStorage.^2+modelLoss.^2));
                                    modelErrorTime = sse_global(visco.forces,...
                                        LR_Maxwell(tempResults.bestParams{k},visco.times,visco.dts,visco.indentations,visco.tipSize,visco.nu,visco.tipGeom,tempResults.elasticSetting,tempResults.fluidSetting));
                                case 'mse'
                                    modelErrorHarmonic = mse_global(simAbsMod,sqrt(modelStorage.^2+modelLoss.^2),numel(harmonicSettings.bestParams));
                                    modelErrorTime = mse_global(visco.forces,...
                                        LR_Maxwell(tempResults.bestParams{k},visco.times,visco.dts,visco.indentations,visco.tipSize,visco.nu,visco.tipGeom,tempResults.elasticSetting,tempResults.fluidSetting),...
                                        numel(harmonicSettings.bestParams)*length(tipSize));
                            end
                            
                        case 'voigt'
                            recordColNum = 2;
                            [modelStorage,modelLoss,~] = visco.harmonics_Voigt(omega,harmonicSettings);
                            switch errorfunc
                                case 'sse'
                                    modelErrorHarmonic = sse_global(simAbsMod,sqrt(modelStorage.^2+modelLoss.^2));
                                    modelErrorTime = sse_global(visco.indentations,...
                                        LR_Voigt(tempResults.bestParams{k},visco.times,visco.dts,visco.forces,visco.tipSize,visco.nu,visco.tipGeom,tempResults.elasticSetting,tempResults.fluidSetting));
                                case 'mse'
                                    modelErrorHarmonic = mse_global(simAbsMod,sqrt(modelStorage.^2+modelLoss.^2),numel(harmonicSettings.bestParams));
                                    modelErrorTime = mse_global(visco.indentations,...
                                        LR_Voigt(tempResults.bestParams{k},visco.times,visco.dts,visco.forces,visco.tipSize,visco.nu,visco.tipGeom,tempResults.elasticSetting,tempResults.fluidSetting),...
                                        numel(harmonicSettings.bestParams)*length(tipSize));
                            end
                            
                        case 'plr'
                            recordColNum = 3;
                            harmonicSettings.dt = dt;
                            harmonicSettings.nu_sample = mode(cell2mat(cellfun(@(x)mode(round(x,4,'significant')),nu,'UniformOutput',false)));
                            [modelStorage,modelLoss,~] = visco.harmonics_PLR(omega,harmonicSettings);
                            switch errorfunc
                                case 'sse'
                                    modelErrorHarmonic = sse_global(simAbsMod,sqrt(modelStorage.^2+modelLoss.^2));
                                    modelErrorTime = sse_global(visco.indentations,...
                                        LR_PLR(tempResults.bestParams{1},visco.times,visco.dts,visco.indentations,visco.tipSize,visco.nu,visco.tipGeom,tempResults.elasticSetting,tempResults.fluidSetting));
                                case 'mse'
                                    modelErrorHarmonic = mse_global(simAbsMod,sqrt(modelStorage.^2+modelLoss.^2),numel(harmonicSettings.bestParams));
                                    modelErrorTime = mse_global(visco.indentations,...
                                        LR_PLR(tempResults.bestParams{1},visco.times,visco.dts,visco.indentations,visco.tipSize,visco.nu,visco.tipGeom,tempResults.elasticSetting,tempResults.fluidSetting),...
                                        numel(harmonicSettings.bestParams)*length(tipSize));
                            end
                            
                        otherwise
                            error('The model in your results structure was not recognized.')
                    end
                    
                    % Decide whether to save the results for the plot
                    switch errortype
                        case 'time'
                            if loopErrorRecord(1,recordColNum) > modelErrorTime
                                barDataLoop(1,recordColNum) = log10(modelErrorTime);
                                numArmLoop(1,recordColNum) = k;
                                if k == simSettings.n_elements
                                    armLabelColorLoop{1,recordColNum} = 'g';
                                else
                                    armLabelColorLoop{1,recordColNum} = 'r';
                                end
                                loopErrorRecord(1,recordColNum) = modelErrorTime;
                            end
                        case 'harmonic'
                            if loopErrorRecord(1,recordColNum) > modelErrorHarmonic
                                barDataLoop(1,recordColNum) = log10(modelErrorHarmonic);
                                numArmLoop(1,recordColNum) = k;
                                if k == simSettings.n_elements
                                    armLabelColorLoop{1,recordColNum} = 'g';
                                else
                                    armLabelColorLoop{1,recordColNum} = 'r';
                                end
                                loopErrorRecord(1,recordColNum) = modelErrorHarmonic;
                            end
                    end
                    
                end

            end

        end
        
        % Stack results onto the running list
        barData = vertcat(barData,barDataLoop);
        numArmData = vertcat(numArmData,numArmLoop);
        armLabelColorData = vertcat(armLabelColorData,armLabelColorLoop);
        testCondList = vertcat(testCondList,simNum);
        testCondData = vertcat(testCondData,{int2str(simNum)});

    end
    
    % Make the bar chart for the current subplot
    hold on
    [~,sortidx] = sort(testCondList);
    X = categorical(testCondData(sortidx));
    X = reordercats(X,testCondData(sortidx));
    b = bar(X,barData(sortidx,:));
    for i = 1:size(barData,2)
%         b(i).LineStyle = linestyleList{i};
        xtip = b(i).XEndPoints;
        ytip = b(i).YEndPoints;
        barLabel = string(numArmData(:,i))';
        textColor = armLabelColorData{:,i};
        htext = text(xtip,ytip,barLabel,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'FontWeight','bold');
        for j = 1:length(htext)
            if strcmp(armLabelColorData{j,i},'r')
                htext(j).Color = [0.6350 0.0780 0.1840];
            elseif strcmp(armLabelColorData{j,i},'g')
                htext(j).Color = [0.4660 0.6740 0.1880];
            end
        end
    end
    hold off

end

% Decide whether to save the results for the plot
switch errortype
    case 'time'
        switch errorfunc
            case 'sse'
                ylabelStr = sprintf('Time Series Error log10(SSE)');
            case 'mse'
                ylabelStr = sprintf('Time Series Error log10(MSE)');
        end
    case 'harmonic'
        switch errorfunc
            case 'sse'
                ylabelStr = sprintf('$\\tilde{E}$ Error log10(SSE)');
            case 'mse'
                ylabelStr = sprintf('$\\tilde{E}$ Error log10(MSE)');
        end
end

switch lower(solverChoice)
    case 'nelder-mead'
        savePrep = 'Nelder';
    case 'annealing'
        savePrep = 'Annealing';
    case 'nls'
        savePrep = 'NLS';
end

figure(errorPlot)
Nplots = numel(findall(errorPlot.Children,'type','axes'));
if Nplots > 1
    nRows = ceil(Nplots/nCols);
    for i = 1:Nplots
        subplot(nRows,nCols,i)
        tempTitle = strsplit(pathList{i},{'/','\'});
        title(sprintf('%s',tempTitle{end}))
        if any(i==(1:nCols:Nplots))
            ylabel(ylabelStr,'fontweight','bold','interpreter','latex')
        end
        if any(i>(nRows-1)*nCols)
            xlabel('Test Condition No.','fontweight','bold','interpreter','latex')
        end
        legend({'Gen. Maxwell', 'Gen. Voigt', 'PLR'},'location','best')
        grid off
        box on
        
        % Add margin to Y scaling
        limSet = findobj(gca, '-property', 'ydata');
        limSet = get(limSet, 'YData');
        limSet = [limSet{:}];
        set(gca,'ylim',[floor(min(limSet))-5 ceil(max(limSet))+5]);
        
    end
    
    sgtitle(sprintf('%s Performance Results',savePrep))
    
else
    title(sprintf('%s Performance Results',savePrep))
    xlabel('Test Condition No.','fontweight','bold')
    ylabel(ylabelStr,'fontweight','bold')
    legend({'Gen. Maxwell', 'Gen. Voigt', 'PLR'},'location','best')
    grid off
    box on
    
    % Add margin to Y scaling
    limSet = findobj(gca, '-property', 'ydata');
    limSet = get(limSet, 'YData');
    limSet = [limSet{:}];
    set(gca,'ylim',[floor(min(limSet))-5 ceil(max(limSet))+5]);
end

figure(errorPlot)
saveas(errorPlot,[outputpathname sprintf('\\%s-ErrorResults-BarComparison.fig',savePrep)])
% saveas(errorPlot,[outputpathname sprintf('\\%s-ErrorResults-BarComparison.jpg',savePrep)])
print(errorPlot,[outputpathname sprintf('\\%s-ErrorResults-BarComparison.png',savePrep)],'-dpng','-r600');
    
% Open the originally requested directory
winopen(outputpathname);