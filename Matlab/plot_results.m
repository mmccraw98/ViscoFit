clear all
close all
clc

addpath(genpath([pwd filesep 'lib']));
addpath(genpath([pwd filesep 'plotting']));

stillRunning = true;
while stillRunning
    
    % User-Defined Settings
    nCols = 2;
    errortype = 'mse';
    figX = 400;
    figY = 150;
    figWid = 800;
    figHeight = 800;
    
    if exist('originalPath','var')
        startdir = originalPath;
    else
        startdir = pwd;
    end
    
    % Pick the AFM Data Directory and Choose Data Extraction Settings
    originalPath = uigetdir(startdir,...
            'Select the Folder Containing Your AFM Files');

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

    % Clear old figures if they exist
    if ~exist('timeErrorPlotNelder','var')
        timeErrorPlotNelder = figure('Position',[figX figY figWid figHeight]);
    else
        try
            figure(timeErrorPlotNelder)
            clf
        catch
            clearvars timeErrorPlotNelder
            timeErrorPlotNelder = figure('Position',[figX figY figWid figHeight]);
        end
    end

    if ~exist('timeErrorPlotAnnealing','var')
        timeErrorPlotAnnealing = figure('Position',[figX+100 figY figWid figHeight]);
    else
        try
            figure(timeErrorPlotAnnealing)
            clf
        catch
            clearvars timeErrorPlotAnnealing
            timeErrorPlotAnnealing = figure('Position',[figX+100 figY figWid figHeight]);
        end
    end

    if ~exist('timeErrorPlotNLS','var')
        timeErrorPlotNLS = figure('Position',[figX+200 figY figWid figHeight]);
    else
        try
            figure(timeErrorPlotNLS)
            clf
        catch
            clearvars timeErrorPlotNLS
            timeErrorPlotNLS = figure('Position',[figX+200 figY figWid figHeight]);
        end
    end

    % Make a list of unique colors for the test conditions
    testConditionColors = num2cell(jet(numel(Folders)),2);

    % Pre-process some results to detect how many subplots we need
    n_elements_tot = [];
    for i = 1:length(Folders)
        path = Folders{i};
        Files = dir([path '\FitResults*.mat']);
        for j = 1:length(Files)
            tempStr = cell2mat(regexp(Files(j).folder,'TestCondition[0-9]*','match'));
            simNum_temp = str2double(cell2mat(regexp(tempStr,'\d*','match')));
            simParams_temp = getSimParams(simNum_temp);
            simParams_temp = horzcat(simParams_temp(1),0,simParams_temp(2:end));
            n_elements_tot = horzcat(n_elements_tot,(numel(simParams_temp)-2)/2);
        end
    end
    n_subplots = numel(unique(n_elements_tot));
    nRows = ceil(n_subplots/nCols);
    
    % Define error functions
    sse_global = @(data,model) sum((data-model).^2,'all');
    mse_global = @(data,model,n) sum(((data-model).^2)./(movvar(data,3).^2),'all')./(length(data)-n);

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

        for j_dir = 1:length(Files)
            fitResults = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
            dataLabels = fieldnames(fitResults);
            resultsContainer = cell(numel(dataLabels),1);
            for j = 1:numel(dataLabels)
                resultsContainer{j} = fitResults.(dataLabels{j});
            end

            % Get simulation parameters
            tempStr = cell2mat(regexp(Files(j_dir).folder,'TestCondition[0-9]*','match'));
            simNum = str2double(cell2mat(regexp(tempStr,'\d*','match')));
            simParams = getSimParams(simNum);
            simParams = horzcat(simParams(1),0,simParams(2:end));

            simSettings = struct;
            simSettings.model = 'maxwell';
            simSettings.bestParams = simParams;
            simSettings.elasticSetting = true;
            simSettings.fluidSetting = false;
            simSettings.n_elements = (numel(simParams)-2)/2;

            markerStyles = {'o','d','s'};
            colorStyles = {'r','b','g','m'};

            for j = 1:numel(resultsContainer)
                tempResults = resultsContainer{j};

                tempSettings = struct;
                tempSettings.minTimescale = 1e-4;
                tempSettings.nu = nu;
                tempSettings.tipGeom = "spherical";
                tempSettings.fitLog = false;
                visco = ViscoFit(forces,times,indentations,tipSize,tempSettings);
                [simStorage,simLoss,simAng] = visco.harmonics_Maxwell(omega,simSettings);
                simAbsMod = sqrt(simStorage.^2+simLoss.^2);

                for k = 1:numel(tempResults.bestParams)

                    harmonicSettings = struct;
                    harmonicSettings.bestParams = tempResults.bestParams{k};
                    harmonicSettings.elasticSetting = tempResults.elasticSetting;
                    harmonicSettings.fluidSetting = tempResults.fluidSetting;
                    harmonicSettings.model = tempResults.model;

                    markerSize = 100+numel(tempResults.bestParams)*100*(1-k/numel(tempResults.bestParams));

                    switch lower(tempResults.solver)
                        case 'nelder-mead'
                            figure(timeErrorPlotNelder)
                            if n_subplots > 1
                                subplot(nRows,nCols,simSettings.n_elements)
                            end
                        case 'annealing'
                            figure(timeErrorPlotAnnealing)
                            if n_subplots > 1
                                subplot(nRows,nCols,simSettings.n_elements)
                            end
                        case 'nls'
                            figure(timeErrorPlotNLS)
                            if n_subplots > 1
                                subplot(nRows,nCols,simSettings.n_elements)
                            end
                        otherwise
                            error('The solver in your results structure was not recognized.')
                    end

                    switch lower(tempResults.model)
                        case 'maxwell'
                            markerType = markerStyles{1};
                            [modelStorage,modelLoss,~] = visco.harmonics_Maxwell(omega,harmonicSettings);
                            switch errortype
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
                            markerType = markerStyles{2};
                            [modelStorage,modelLoss,~] = visco.harmonics_Voigt(omega,harmonicSettings);
                            switch errortype
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
                            markerType = markerStyles{3};
                            harmonicSettings.dt = dt;
                            harmonicSettings.nu_sample = mode(cell2mat(cellfun(@(x)mode(round(x,4,'significant')),nu,'UniformOutput',false)));
                            [modelStorage,modelLoss,~] = visco.harmonics_PLR(omega,harmonicSettings);
                            switch errortype
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

                    hold on
                    scatter(modelErrorTime,modelErrorHarmonic,markerSize,markerType,'MarkerFaceColor',...
                        testConditionColors{i_dir},'MarkerFaceAlpha',1,'MarkerEdgeColor','k')
                    hold off

                end

            end

        end

    end

    switch errortype
        case 'sse'
            ylabelStr = sprintf('Absolute Modulus Error (SSE)');
            xlabelStr = sprintf('Time Series Error (SSE)');
        case 'mse'
            ylabelStr = sprintf('Absolute Modulus Error (MSE)');
            xlabelStr = sprintf('Time Series Error (MSE)');
    end
    
    figure(timeErrorPlotNelder)
    Nplots = numel(findall(timeErrorPlotNelder.Children,'type','axes'));
    if Nplots > 1
        nRows = ceil(Nplots/nCols);
        for i = 1:Nplots
            subplot(nRows,nCols,i)
            title(sprintf('%d-Term Samples',i))
            if any(i==(1:nCols:Nplots))
                ylabel(ylabelStr,'fontweight','bold')
            end
            if any(i>(nRows-1)*nCols)
                xlabel(xlabelStr,'fontweight','bold')
            end
            set(gca,'XScale','log','YScale','log')
            if strcmp(findobj(gca, '-property', 'xscale'),'log')
                limSet = findobj(gca, '-property', 'xdata');
                limSet = get(limSet, 'XData');
                limSet = [limSet{:}];
                set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
            end
            grid on
        end
        sgtitle('Nelder-Mead Performance Results')
    else
        title('Nelder-Mead Performance Results')
        xlabel(xlabelStr,'fontweight','bold')
        ylabel(ylabelStr,'fontweight','bold')
        set(gca,'XScale','log','YScale','log')
        if strcmp(findobj(gca, '-property', 'xscale'),'log')
            limSet = findobj(gca, '-property', 'xdata');
            limSet = get(limSet, 'XData');
            limSet = [limSet{:}];
            set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
        end
        grid on
    end

    figure(timeErrorPlotAnnealing)
    Nplots = numel(findall(timeErrorPlotAnnealing.Children,'type','axes'));
    if Nplots > 1
        nRows = ceil(Nplots/nCols);
        for i = 1:Nplots
            subplot(nRows,nCols,i)
            title(sprintf('%d-Term Samples',i))
            if any(i==(1:nCols:Nplots))
                ylabel(ylabelStr,'fontweight','bold')
            end
            if any(i>(nRows-1)*nCols)
                xlabel(xlabelStr,'fontweight','bold')
            end
            set(gca,'XScale','log','YScale','log')
            if strcmp(findobj(gca, '-property', 'xscale'),'log')
                limSet = findobj(gca, '-property', 'xdata');
                limSet = get(limSet, 'XData');
                limSet = [limSet{:}];
                set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
            end
            grid on
        end
        sgtitle('Simulated Annealing Performance Results')
    else
        title('Simulated Annealing Performance Results')
        xlabel(xlabelStr,'fontweight','bold')
        ylabel(ylabelStr,'fontweight','bold')
        set(gca,'XScale','log','YScale','log')
        if strcmp(findobj(gca, '-property', 'xscale'),'log')
            limSet = findobj(gca, '-property', 'xdata');
            limSet = get(limSet, 'XData');
            limSet = [limSet{:}];
            set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
        end
        grid on
    end

    figure(timeErrorPlotNLS)
    Nplots = numel(findall(timeErrorPlotNLS.Children,'type','axes'));
    if Nplots > 1
        nRows = ceil(Nplots/nCols);
        for i = 1:Nplots
            subplot(nRows,nCols,i)
            title(sprintf('%d-Term Samples',i))
            if any(i==(1:nCols:Nplots))
                ylabel(ylabelStr,'fontweight','bold')
            end
            if any(i>(nRows-1)*nCols)
                xlabel(xlabelStr,'fontweight','bold')
            end
            set(gca,'XScale','log','YScale','log')
            if strcmp(findobj(gca, '-property', 'xscale'),'log')
                limSet = findobj(gca, '-property', 'xdata');
                limSet = get(limSet, 'XData');
                limSet = [limSet{:}];
                set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
            end
            grid on
        end
        sgtitle('Nonlinear Least-Squares Performance')
    else
        title('Nonlinear Least-Squares Performance')
        xlabel(xlabelStr,'fontweight','bold')
        ylabel(ylabelStr,'fontweight','bold')
        set(gca,'XScale','log','YScale','log')
        if strcmp(findobj(gca, '-property', 'xscale'),'log')
            limSet = findobj(gca, '-property', 'xdata');
            limSet = get(limSet, 'XData');
            limSet = [limSet{:}];
            set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
        end
        grid on
    end

    figure(timeErrorPlotNelder)
    saveas(timeErrorPlotNelder,[originalPath '\ErrorResults-NelderMead.fig'])
%     saveas(timeErrorPlotNelder,[originalPath '\ErrorResults-NelderMead.jpg'])
    print(timeErrorPlotNelder,[originalPath '\ErrorResults-NelderMead.png'],'-dpng','-r300');
    figure(timeErrorPlotAnnealing)
    saveas(timeErrorPlotAnnealing,[originalPath '\ErrorResults-Annealing.fig'])
%     saveas(timeErrorPlotAnnealing,[originalPath '\ErrorResults-Annealing.jpg'])
    print(timeErrorPlotAnnealing,[originalPath '\ErrorResults-Annealing.png'],'-dpng','-r300');
    figure(timeErrorPlotNLS)
    saveas(timeErrorPlotNLS,[originalPath '\ErrorResults-NLS.fig'])
%     saveas(timeErrorPlotNLS,[originalPath '\ErrorResults-NLS.jpg'])
    print(timeErrorPlotNLS,[originalPath '\ErrorResults-NLS.png'],'-dpng','-r300');

    % Prompt user
    answer = questdlg('Would you like to analyze another directory?', ...
        'Options', ...
        'No','Yes','Yes');

    % Handle response
    switch answer
        case 'No'
            clearvars -except originalPath
            stillRunning = false;
        case 'Yes'
            close all
            clearvars -except stillRunning originalPath
    end
    
end
    
% Open the originally requested directory
winopen(originalPath);