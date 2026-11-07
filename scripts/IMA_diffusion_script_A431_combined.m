% This script is used to estimate diffusion constants from TR structures.

%% Define some directories and parameters.
FrameRate = 18.18;        % frames per second
PixelSize = 0.0677;   % um % 1X83, with 60X oil and 1.6X mag. 
% PixelSize = 0.1667; % ix71
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FrameRate = FrameRate;
SMF.Data.PixelSize = PixelSize;
% DMinMax = [0.001, 2];
DMinMax = [0, inf]; % um^2/s, [min., max.] for MSD results (cell and trajectory)
NComponents = 3; % number of diffusing components for CDF/MLE fitting
% NComponents = 2; % number of diffusing components for CDF/MLE fitting
% Minimum and maximum frame numbers for analysis.
% NOTE: if MinFrameNum < 0 (e.g., -3), then it is the absolute number of
%       frames (e.g., +3) to retain per trajectory.
%       if MaxFrameNum > 0 and < 1 (e.g., 0.3), then it is the fraction of the
%       frames (e.g., 30%) to retain per trajectory.
% MinFrameNum = 1;   MaxFrameNum = Inf;   % all frames from 1 to Inf per traj
%MinFrameNum = -10;   MaxFrameNum = 50;  % first 10 frames per trajectory
%MinFrameNum = 1;   MaxFrameNum = 0.2;   % first 20% of the frames per traj
% Specify the maximum starting frame number for a trajectory to be included in
% the diffusion analyses.
% MaxStartingFrameNum = 100;
%MaxStartingFrameNum = 3;
% 

BaseDir = 'O:\Cell Path\Lidke Lab\Genmab\IMAdams\DATA\Combined_A431_0904-0915';
experimentName =  'A431_IgG-488_dark-IgG_C1Q_n3_appended'; % 'Qdot655';
chName = '1'; % Qdot-655 -> '2'; Qdot-585 -> '1'
plotTitle = 'A431 n=3 experiments (old C1q)';

ConditionNames = {...
    '2F8'...
    '+dark 2F8'...
    '2F8 + C1q'...
    'E430G'...
    '+dark E430G'...
    'E430G + C1q'...
    'E345R'...
    '+dark E345R'...
    'E345R + C1q'...
    'RGY'...
    '+dark RGY'...
    'RGY + C1q'...
    };

NConditions = numel(ConditionNames);

SaveDir = fullfile(BaseDir,strcat('D-analysis_', experimentName, '_ch', chName));
% SaveDir = 'O:\Cell Path\Lidke Lab\IMAdams\Data\CHO-HA-EGFR-L858R\Single particle tracking\20251016_CHO_WThv2_serumpilot\Results1020\D_10202025';

filterBy = strcat(chName,'_Results.mat');

% filterBy = '2_Results.mat';

% % load file list struct
% FileList = load('O:\Cell Path\Lidke Lab\Genmab\IMAdams\DATA\Combined_A431_0904-0915\D-analysis_A431_IgG-488_dark-IgG_C1Q_n3_ch1\A431_IgG-488_dark-IgG_C1Q_n3FileList.mat');
% FileList = FileList.FileList; % struct -> cell


% Run the file picker
FileList = cell(NConditions, 1);
for ii = 1:NConditions
    fprintf('Pick files for condition %s\n', ConditionNames{ii})
    FileList{ii} = uipickfiles('FilterSpec', BaseDir,...
        'REFilter', filterBy, ...
        'Prompt', ...
        sprintf('Pick files for condition %s\n', ConditionNames{ii}));
end

% Quality control: check for duplicate files and report files per group.
FileListReport = checkFileList(FileList, ConditionNames);

% If SaveDir is non-empty but doesn't exist yet, create the directory.
if ~isempty(SaveDir) && ~isfolder(SaveDir)
    mkdir(SaveDir);
end
fileString = fullfile(SaveDir, strcat(experimentName, "FileList"));
save(fileString, "FileList");

%% Loop through the TR structures and estimate the diffusion constants.
DEstimator = smi_stat.DiffusionEstimator([], SMF);
DEstimator.UnitFlag = true;
DEstimator.FrameLagRange = [1, 20];   % Range of frame lags used to estimate D (Default = [1, 5])
DEstimator.NFitPoints = 5; % Number of MSD points to be fit (scalar, integer)(Default = 5)
DEstimator.FitMethod = 'WeightedLS';
DEstimator.FitTarget = 'MSD';
DEstimator.Verbose = 1; 
% Indexed by {condition #}
MSDEnsemble = cell(NConditions, 1);
DAllTraj = cell(NConditions, 1);
DAllCells = cell(NConditions, 1);
MeanPhotonsAllTraj = cell(NConditions, 1);
DEnsemble = zeros(NConditions, 1);
DEnsembleSE = DEnsemble;
KeptFiles = cell(NConditions, 1);
NCells = zeros(NConditions, 1);
NTraj = zeros(NConditions, 1);
NTrajEnsemble = zeros(NConditions, 1);
NJumps = zeros(NConditions, 1);
NJumpsFit = zeros(NConditions, 1);
for dd = 1:NConditions
    fprintf('Analyzing condition %i of %i...\n', dd, NConditions)
    
    % Load the TR structures for this condition.
    NTRCurrent = numel(FileList{dd});
    TRAll = struct([]);
    for ii = 1:NTRCurrent
        % Load the TR structure and modify some properties as appropriate.
        load(FileList{dd}{ii}, 'TR')
        if ~isfield(TR, 'ConnectID')
            for jj = 1:numel(TR)
                TR(jj).ConnectID = TR(jj).TrajectoryID;
            end
        end

        % % Choose a frame number window from the TR structure. Specify
        % MinFrameNum, MaxFrameNum
        
        % TR = smi_core.TrackingResults.windowTR(TR, MinFrameNum, MaxFrameNum);
        % TR = smi_core.TrackingResults.windowStartTR(TR, MaxStartingFrameNum);
        % 
        % Estimate the trajectory-wise diffusion constants (one diffusion
        % constant per trajectory, excluding trajectories that were too
        % short).

        % skip analyzing empty results 
        if (numel(TR)<2)
            continue
        end

        DEstimator.TR = TR;
        DiffusionStruct = DEstimator.estimateDiffusionConstant();

        % Determine if this cell should be excluded based on our
        % thresholds.
        if ((DiffusionStruct(2).DiffusionConstant<DMinMax(1)) ...
                || (DiffusionStruct(2).DiffusionConstant>DMinMax(2)))
            % Don't store the results from this cell.
            continue
        end
        
        % Store the current results in the concatenated arrays.
        KeptFiles{dd} = [KeptFiles{dd}; FileList{dd}(ii)];    
        DAllCells{dd} = [DAllCells{dd}; ...
            [DiffusionStruct(2).DiffusionConstant, ...
            DiffusionStruct(2).DiffusionConstantSE]];
        TRAll = smi_core.TrackingResults.catTR(TRAll, TR);
    end
    
    % Find diffusion constants for the entire condition.
    DEstimator.TR = TRAll;
    DiffusionStruct = DEstimator.estimateDiffusionConstant();
    MSDEnsemble{dd} = DEstimator.MSDEnsemble;
    DAllTrajNew = [DiffusionStruct(1).DiffusionConstant, ...
        DiffusionStruct(1).DiffusionConstantSE];
    KeepBool = ((DAllTrajNew(:, 1)>=DMinMax(1)) ...
        & (DAllTrajNew(:, 1)<=DMinMax(2)));
    DAllTraj{dd} = [DAllTraj{dd}; DAllTrajNew(KeepBool, :)];
    NTraj(dd) = sum(KeepBool);
    % Compute average photon counts per trajectory.
    MeanPhotonsAllTrajNew = ...
        arrayfun(@(i) mean(TRAll(i).Photons), 1:numel(TRAll));
    MeanPhotonsAllTraj{dd} = MeanPhotonsAllTrajNew(KeepBool);
    DEnsemble(dd) = DiffusionStruct(2).DiffusionConstant;
    DEnsembleSE(dd) = DiffusionStruct(2).DiffusionConstantSE;
    NCells(dd) = numel(KeptFiles{dd});
    NTrajEnsemble(dd) = numel(TRAll);
    NJumps(dd) = sum(DEstimator.MSDEnsemble.NPoints);
    NJumpsFit(dd) = sum(DEstimator.MSDEnsemble.NPoints(DEstimator.NFitPoints));
end

% Repeat above but estimate D from CDF of jumps (instead of fitting an MSD)
DEstimator = smi_stat.DiffusionEstimator([], SMF);
DEstimator.UnitFlag = true;
DEstimator.FrameLagRange = [2,2];
DEstimator.NFitPoints = inf;
DEstimator.FitMethod = 'LS';
DEstimator.FitTarget = 'CDFOfJumps';
DEstimator.DiffusionModel = 'Brownian';
DEstimator.NComponents = NComponents;
DEstimator.FitIndividualTrajectories = false;
DEstimator.EstimateSEs = false;
CDFEnsemble = cell(NConditions, 1);
DiffusionStruct = cell(NConditions, 1);
for dd = 1:NConditions
    fprintf('Analyzing condition %i of %i...\n', dd, NConditions)
    
    % Load the TR structures for this condition.
    NTRCurrent = numel(FileList{dd});
    TRAll = struct([]);
    for ii = 1:NTRCurrent
        % Load the TR structure and modify some properties as appropriate.
        load(FileList{dd}{ii}, 'TR')
        if ~isfield(TR, 'ConnectID')
            for jj = 1:numel(TR)
                TR(jj).ConnectID = TR(jj).TrajectoryID;
            end
        end

        % Choose a frame number window from the TR structure.
        % TR = smi_core.TrackingResults.windowTR(TR, MinFrameNum, MaxFrameNum);
        % TR = smi_core.TrackingResults.windowStartTR(TR, MaxStartingFrameNum);
        % 
        % Concatenate the TR into TRAll.
        TRAll = smi_core.TrackingResults.catTR(TRAll, TR);
    end
    
    % Find diffusion constants for the entire condition.
    DEstimator.TR = TRAll;
    DiffusionStruct{dd} = DEstimator.estimateDiffusionConstant();
    CDFEnsemble{dd} = DEstimator.MSDEnsemble;
    NTrajEnsemble(dd) = numel(TRAll);
    NJumps(dd) = numel(DEstimator.MSDEnsemble.CDFOfJumps);
end
%% Save key variables for later reload.
AnalysisFile = fullfile(SaveDir, 'DiffusionAnalysisData.mat');
SaveVars = {'FileList', 'KeptFiles', 'ConditionNames', 'NConditions', ...
    'BaseDir', 'experimentName', 'chName', 'plotTitle', 'SaveDir', ...
    'FrameRate', 'PixelSize', 'DMinMax', 'NComponents', 'SMF', ...
    'MSDEnsemble', ...
    'DAllTraj', 'DAllCells', 'MeanPhotonsAllTraj', ...
    'DEnsemble', 'DEnsembleSE', ...
    'NCells', 'NTraj', 'NTrajEnsemble', 'NJumps', 'NJumpsFit', ...
    'CDFEnsemble', 'DiffusionStruct'};
% Save as compressed -v7 (smaller/faster); fall back to -v7.3 (HDF5) if any
% variable exceeds the 2 GB limit of -v7.
lastwarn('');
save(AnalysisFile, SaveVars{:});
if ~isempty(lastwarn)
    fprintf('-v7 save warned (%s); re-saving as -v7.3.\n', lastwarn);
    save(AnalysisFile, SaveVars{:}, '-v7.3');
end
fprintf('Saved diffusion analysis data to %s\n', AnalysisFile);
%% Plot some results.

% Plot the MSD for all conditions (ensemble over trajectories and cells).
% PlotAxes = axes(figure());
% hold(PlotAxes, 'on')
% LegendEntries = cell(NConditions, 1);
% ConditionColors = lines(NConditions);
% for nn = 1:NConditions
%     plot(PlotAxes, MSDEnsemble{nn}.FrameLags/SMF.Data.FrameRate, ...
%         MSDEnsemble{nn}.MSD*(SMF.Data.PixelSize^2), ...
%         'LineWidth', 2, 'Color', ConditionColors(nn, :))
%     LegendEntries{nn} = sprintf('%s', ...
%         ConditionNames{nn}); % remove D from legend
%     LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
%         ConditionNames{nn}, DEnsemble(nn), 1.96*DEnsembleSE(nn));
% end
% title(PlotAxes, strcat(plotTitle, ' MSD'))
% PlotAxes.FontWeight = 'bold';
% PlotAxes.XLim = [0, 1];
% legend(PlotAxes, LegendEntries, 'Location', 'northwest')
% xlabel(PlotAxes, '\Deltat (s)', 'Interpreter', 'tex', 'FontWeight', 'bold')
% ylabel(PlotAxes, 'MSD (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
% if ~isempty(SaveDir)
%    saveas(gcf, fullfile(SaveDir, 'MSDvsDeltat.png'));
%    saveas(gcf, fullfile(SaveDir, 'MSDvsDeltat.fig'));
% end
%%
NGroups = 4;
NSubConditions = 3;

PlotAxes = axes(figure());
hold(PlotAxes, 'on')
LegendEntries = cell(NConditions, 1);

GroupColors = lines(NGroups);              % one color per group
LineStyles = {'-', '--', ':'};             % one linestyle per sub-condition

for nn = 1:NConditions
    GroupIndex = ceil(nn / NSubConditions);
    SubIndex = mod(nn - 1, NSubConditions) + 1;

    plot(PlotAxes, MSDEnsemble{nn}.FrameLags/SMF.Data.FrameRate, ...
        MSDEnsemble{nn}.MSD*(SMF.Data.PixelSize^2), ...
        'LineWidth', 2, ...
        'Color', GroupColors(GroupIndex, :), ...
        'LineStyle', LineStyles{SubIndex})

    LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
        ConditionNames{nn}, DEnsemble(nn), 1.96*DEnsembleSE(nn));
end

title(PlotAxes, strcat(plotTitle, ' MSD'))
PlotAxes.FontWeight = 'bold';
PlotAxes.XLim = [0, 1];
legend(PlotAxes, LegendEntries, 'Location', 'northwest')
xlabel(PlotAxes, '\Deltat (s)', 'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxes, 'MSD (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'MSDvsDeltat.png'));
   saveas(gcf, fullfile(SaveDir, 'MSDvsDeltat.fig'));
end
%%
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
LegendEntries = cell(NConditions, 1);
ConditionColors = lines(NConditions);
for nn = 1:NConditions
    tLags = MSDEnsemble{nn}.FrameLags / SMF.Data.FrameRate;
    MSDPlot = MSDEnsemble{nn}.MSD * (SMF.Data.PixelSize^2);

    % Compute SEM at each lag from raw squared displacements.
    NLags = numel(MSDEnsemble{nn}.FrameLags);
    SEM = zeros(NLags, 1);
    for ii = 1:NLags
        SD = MSDEnsemble{nn}.SquaredDisplacement( ...
            MSDEnsemble{nn}.FrameLagsAll == MSDEnsemble{nn}.FrameLags(ii));
        SEM(ii) = std(SD) / sqrt(numel(SD));
    end
    SEMPlot = SEM * (SMF.Data.PixelSize^2);

    % Shading behind t
    XPatch = [tLags; flip(tLags)];
    YPatch = [MSDPlot + SEMPlot; flip(MSDPlot - SEMPlot)];              
    patch(PlotAxes, XPatch, YPatch, ConditionColors(nn, :), ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    plot(PlotAxes, tLags, MSDPlot, ...
        'LineWidth', 2, 'Color', ConditionColors(nn, :))
    LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
         ConditionNames{nn}, DEnsemble(nn), 1.96*DEnsembleSE(nn));
end
PlotAxes.FontWeight = 'bold';                                       
PlotAxes.XLim = [0, 1]
legend(PlotAxes, LegendEntries, 'Location', 'northwest')
xlabel(PlotAxes, '\Deltat (s)', 'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxes, 'MSD (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'MSDvsDeltatSHADED.png'));
   saveas(gcf, fullfile(SaveDir, 'MSDvsDeltatSHADED.fig'));
end


% Plot the CDFs for all conditions (ensemble over trajectories and cells).
SqDisp = cell(NConditions, 1);
CDFs = cell(NConditions, 1);
for nn = 1:NConditions
    SqDisp{nn} = CDFEnsemble{nn}.SortedSquaredDisp * (SMF.Data.PixelSize^2);
    CDFs{nn} = CDFEnsemble{nn}.CDFOfJumps;
end
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
ConditionColors = lines(NConditions);
smi_vis.plotCDFs(PlotAxes, SqDisp, CDFs, ConditionColors)
PlotAxes.FontWeight = 'bold';
axis(PlotAxes, 'tight')
PlotAxes.XScale = 'log';
legend(PlotAxes, ConditionNames, 'Location', 'best')
xlabel(PlotAxes, 'r^2 (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxes, 'CDF(r^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
title(PlotAxes, strcat(plotTitle, ' CDF'))
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'CDFvsr2.png'));
   saveas(gcf, fullfile(SaveDir, 'CDFvsr2.fig'));
end

% Plot the CDFs for all conditions including the fits.
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
ModelCDFs = cell(NConditions, 1);
LegendEntries = cell(2*NConditions, 1);
LegendEntries(1:NConditions) = ConditionNames;
for nn = 1:NConditions
    ModelCDFs{nn} = smi_stat.DiffusionEstimator.brownianJumpCDF(...
        [DiffusionStruct{nn}(2).DiffusionConstant, ...
        DiffusionStruct{nn}(2).PopulationRatios], ...
        CDFEnsemble{nn}.SortedSquaredDisp * (SMF.Data.PixelSize^2), ...
        CDFEnsemble{nn}.FrameLags / SMF.Data.FrameRate, ...
        CDFEnsemble{nn}.NPoints, ...
        CDFEnsemble{nn}.LocVarianceSum * (SMF.Data.PixelSize^2));
    LegendEntries{nn + NConditions} = ...
           [sprintf('%s fit', ConditionNames{nn}), newline, ...
            sprintf('D_1 = %.4g, D_2 = %.4g, D_3 = %.4g, Pop_1 = %.4f, Pop_2 = %.4f, Pop_3 = %.4f', ...
               DiffusionStruct{nn}(2).DiffusionConstant, ...
               DiffusionStruct{nn}(2).PopulationRatios) ...
           ];
%   LegendEntries{nn + NConditions} = sprintf('%s fit', ConditionNames{nn});
end
smi_vis.plotCDFs(PlotAxes, SqDisp, CDFs, ConditionColors)
LineHandles = smi_vis.plotCDFs(PlotAxes, SqDisp, ModelCDFs, ConditionColors);
[LineHandles.LineStyle] = deal(':');
PlotAxes.FontWeight = 'bold';
axis(PlotAxes, 'tight')
PlotAxes.XScale = 'log';
legend(PlotAxes, LegendEntries, 'Location', 'best')
xlabel(PlotAxes, 'r^2 (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxes, 'CDF(r^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
title(PlotAxes, sprintf('CDF(r^2) and %i component fits', NComponents))
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'DComponentFits.png'));
   saveas(gcf, fullfile(SaveDir, 'DComponentFits.fig'));
end
%
% Plot the diffusion constants for each cell.
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
DAllCellsIsolated = cellfun(@(X) X(:, 1), DAllCells, ...
    'UniformOutput', false); % isolate D only, not D_SE

% Prepare data for boxplot (needs to be in a specific format)
% Create a grouped vector and corresponding group labels
AllD = vertcat(DAllCellsIsolated{:});
GroupLabels = [];
for nn = 1:NConditions
    GroupLabels = [GroupLabels; repmat(nn, numel(DAllCellsIsolated{nn}), 1)];
end

% Create box plot
boxplot(PlotAxes, AllD, GroupLabels, ...
    'Labels', ConditionNames, ...
    'Colors', 'k', ...
    'Symbol', ''); % suppress outlier markers since we're showing all points

% Overlay individual data points
for nn = 1:NConditions
    % Add jitter to x-positions for visibility
    xJitter = 0.15 * (rand(size(DAllCellsIsolated{nn})) - 0.5);
    scatter(PlotAxes, nn + xJitter, DAllCellsIsolated{nn}, ...
        45, ConditionColors(ceil(nn/3), :), 'filled', 'MarkerFaceAlpha', 0.7)
end
%
% plotSpread(PlotAxes, DAllCellsIsolated, ...
%     'spreadWidth', 0.5, 'distributionColors', ConditionColors)
% for nn = 1:NConditions
%     % Add a horizontal line at the mean value.
%     plot(PlotAxes, nn + 0.25*[-1, 1], mean(DAllCellsIsolated{nn})*[1, 1], ...
%         'Color', [0, 0, 0], 'LineWidth', 2)
% end
PlotAxes.XLim = [0, NConditions+1];
PlotAxes.XTick = 1:NConditions;
PlotAxes.XTickLabel = ConditionNames;
PlotAxes.XTickLabelRotation = 45;
PlotAxes.FontWeight = 'bold';
PlotAxes.YLim = [0,0.16];
PlotAxes.FontSize = 15;
ylabel(PlotAxes, 'D (\mum^2s^{-1})', 'Interpreter', 'tex', ...
    'FontWeight', 'bold')
title(PlotAxes, strcat(plotTitle, ' Cell diffusion coefficients from MSD fits'))
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'DCoeffperCond.png'));
   saveas(gcf, fullfile(SaveDir, 'DCoeffperCond.fig'));
end
%
% Plot the diffusion constants for each trajectory.
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
DAllTrajIsolated = cellfun(@(X) X(:, 1), DAllTraj, ...
    'UniformOutput', false); % isolate D only, not D_SE
plotSpread(PlotAxes, DAllTrajIsolated, ...
    'spreadWidth', 0.5, 'distributionColors', ConditionColors)
for nn = 1:NConditions
    % Add a horizontal line at the mean value.
    plot(PlotAxes, nn + 0.25*[-1, 1], mean(DAllTrajIsolated{nn})*[1, 1], ...
        'Color', [0, 0, 0], 'LineWidth', 2)
end
PlotAxes.XLim = [0, NConditions+1];
PlotAxes.XTick = 1:NConditions;
PlotAxes.XTickLabel = ConditionNames;
PlotAxes.XTickLabelRotation = 45;
PlotAxes.FontWeight = 'bold';
ylabel(PlotAxes, 'D (\mum^2s^{-1})', 'Interpreter', 'tex', ...
    'FontWeight', 'bold')
title(PlotAxes, 'Trajectory diffusion coefficients from MSD fits')
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'TrajDCoeffperCond.png'));
   saveas(gcf, fullfile(SaveDir, 'TrajDCoeffperCond.fig'));
end
%
% Histogram of diffusion coefficients on negative log scale
figure();
PlotAxesHist = axes();
hold(PlotAxesHist, 'on')

% Transform D to -log10(D) for all conditions
DLogAll = cellfun(@(X) -log10(X), DAllTrajIsolated, 'UniformOutput', false);

% Determine bin edges across all data
AllDLog = vertcat(DLogAll{:});
BinEdges = linspace(min(AllDLog), max(AllDLog), 30);

% Plot histograms for each condition
for nn = 1:NConditions
    histogram(PlotAxesHist, DLogAll{nn}, BinEdges, ...
        'FaceColor', ConditionColors(nn, :), ...
        'FaceAlpha', 0.5, ...
        'EdgeColor', 'none', ...
        'Normalization', 'probability', ...
        'DisplayName', ConditionNames{nn})
end

% Format axes
xlabel(PlotAxesHist, '-log_{10}(D) (-log_{10}(\mum^2s^{-1}))', ...
    'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxesHist, 'probability', 'FontWeight', 'bold')
title(PlotAxesHist, 'D from MSD fits ')
legend(PlotAxesHist, 'Location', 'best')
PlotAxesHist.FontWeight = 'bold';

% Add secondary x-axis showing actual D values
ax2 = axes('Position', PlotAxesHist.Position, ...
    'XAxisLocation', 'top', ...
    'Color', 'none', ...
    'YTick', []);
ax2.XLim = PlotAxesHist.XLim;
TickVals = PlotAxesHist.XTick;
ax2.XTick = TickVals;
ax2.XTickLabel = arrayfun(@(x) sprintf('%.3g', 10^(-x)), TickVals, 'UniformOutput', false);
ax2.FontWeight = 'bold';
xlabel(ax2, 'D (\mum^2s^{-1})', 'Interpreter', 'tex', 'FontWeight', 'bold')
if ~isempty(SaveDir)
    saveas(gcf, fullfile(SaveDir, 'TrajDCoeffHist.png'));
    saveas(gcf, fullfile(SaveDir, 'TrajDCoeffHist.fig'));
end

% export N vals
IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, DAllTrajIsolated, fullfile(SaveDir, 'summary.txt'), 'Format', 'detailed')

% Plot the the average photon count for each trajectory.
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
plotSpread(PlotAxes, MeanPhotonsAllTraj, ...
    'spreadWidth', 0.5, 'distributionColors', ConditionColors)
for nn = 1:NConditions
    % Add a horizontal line at the mean value.
    plot(PlotAxes, nn + 0.25*[-1, 1], mean(MeanPhotonsAllTraj{nn})*[1, 1], ...
        'Color', [0, 0, 0], 'LineWidth', 2)
end
PlotAxes.XLim = [0, NConditions+1];
PlotAxes.XTick = 1:NConditions;
PlotAxes.XTickLabel = ConditionNames;
PlotAxes.XTickLabelRotation = 45;
PlotAxes.FontWeight = 'bold';
ylabel(PlotAxes, 'mean photons', 'Interpreter', 'tex', ...
    'FontWeight', 'bold')
title(PlotAxes, 'Trajectory mean photons')
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'TrajPhotonsperCond.png'));
   saveas(gcf, fullfile(SaveDir, 'TrajPhotonsperCond.fig'));
end

% for nn = 1:NConditions
%    figure();
%    hold on
%    plot(MeanPhotonsAllTraj{nn}, DAllTrajIsolated{nn}, 'k.');
%    title(ConditionNames{nn});
%    xlabel('mean photons per trajectory');
%    ylabel('D (\mum^2s^{-1}) per trajectory');
%    hold off
%    if ~isempty(SaveDir)
%       saveas(gcf, fullfile(SaveDir, sprintf('DvsPhotonsCond%d.png', nn)));
%       saveas(gcf, fullfile(SaveDir, sprintf('DvsPhotonsCond%d.fig', nn)));
%    end
% end

%% Save key variables for later reload.
AnalysisFile = fullfile(SaveDir, 'DiffusionAnalysisData.mat');
SaveVars = {'FileList', 'KeptFiles', 'ConditionNames', 'NConditions', ...
    'BaseDir', 'experimentName', 'chName', 'plotTitle', 'SaveDir', ...
    'FrameRate', 'PixelSize', 'DMinMax', 'NComponents', 'SMF', ...
    'MSDEnsemble', ...
    'DAllTraj', 'DAllCells', 'MeanPhotonsAllTraj', ...
    'DEnsemble', 'DEnsembleSE', ...
    'NCells', 'NTraj', 'NTrajEnsemble', 'NJumps', 'NJumpsFit', ...
    'CDFEnsemble', 'DiffusionStruct', ...
    'DFit', 'DFitSE'};
lastwarn('');
save(AnalysisFile, SaveVars{:});
if ~isempty(lastwarn)
    fprintf('-v7 save warned (%s); re-saving as -v7.3.\n', lastwarn);
    save(AnalysisFile, SaveVars{:}, '-v7.3');
end
fprintf('Saved diffusion analysis data to %s\n', AnalysisFile);
%% Reload saved diffusion analysis data.
AnalysisFile = fullfile(SaveDir, 'DiffusionAnalysisData.mat');
load(AnalysisFile)

%% Append new cells to existing groups and update the saved analysis data.
% Run this section on its own (it reloads the saved data itself). For each
% condition you pick the NEW files only (leave the picker empty / cancel to
% skip a condition). New cells are analyzed individually, appended to
% FileList/KeptFiles/DAllCells, and then the condition-level (ensemble, per
% trajectory, and CDF) results are recomputed over old + new cells.
AnalysisFile = fullfile(SaveDir, 'DiffusionAnalysisData.mat');
load(AnalysisFile)
filterBy = strcat(chName, '_Results.mat');

% Back up the existing analysis file before modifying anything.
BackupFile = fullfile(SaveDir, sprintf('DiffusionAnalysisData_backup_%s.mat', ...
    datestr(now, 'yyyymmdd_HHMMSS')));
copyfile(AnalysisFile, BackupFile);
fprintf('Backed up existing analysis data to %s\n', BackupFile);

% Pick the new files for each condition.
NewFileList = cell(NConditions, 1);
for ii = 1:NConditions
    fprintf('Pick NEW files for condition %s\n', ConditionNames{ii})
    Picked = uipickfiles('FilterSpec', BaseDir, 'REFilter', filterBy, ...
        'Prompt', sprintf('Pick NEW files for condition %s', ConditionNames{ii}));
    if isnumeric(Picked) || isempty(Picked)
        Picked = {};   % uipickfiles returns 0 on cancel
    end
    Picked = Picked(:);
    % Drop files that are already part of this condition.
    IsDup = ismember(Picked, FileList{ii});
    if any(IsDup)
        fprintf('  Skipping %d file(s) already in this condition.\n', sum(IsDup))
    end
    NewFileList{ii} = Picked(~IsDup);
end
% Quality control (duplicates across groups / files per group) for new files.
checkFileList(NewFileList, ConditionNames);

% Per-cell and per-condition estimator (same settings as the original run).
DEstimator = smi_stat.DiffusionEstimator([], SMF);
DEstimator.UnitFlag = true;
DEstimator.FrameLagRange = [1, 20];
DEstimator.NFitPoints = 5;
DEstimator.FitMethod = 'WeightedLS';
DEstimator.FitTarget = 'MSD';
DEstimator.Verbose = 1;

% CDF-of-jumps estimator (same settings as the original run).
DEstimatorCDF = smi_stat.DiffusionEstimator([], SMF);
DEstimatorCDF.UnitFlag = true;
DEstimatorCDF.FrameLagRange = [2, 2];
DEstimatorCDF.NFitPoints = inf;
DEstimatorCDF.FitMethod = 'LS';
DEstimatorCDF.FitTarget = 'CDFOfJumps';
DEstimatorCDF.DiffusionModel = 'Brownian';
DEstimatorCDF.NComponents = NComponents;
DEstimatorCDF.FitIndividualTrajectories = false;
DEstimatorCDF.EstimateSEs = false;

for dd = 1:NConditions
    if isempty(NewFileList{dd})
        continue
    end
    fprintf('Appending %d new cell(s) to condition %i of %i (%s)...\n', ...
        numel(NewFileList{dd}), dd, NConditions, ConditionNames{dd})

    % Analyze each new cell individually and append to the group.
    for ii = 1:numel(NewFileList{dd})
        load(NewFileList{dd}{ii}, 'TR')
        if ~isfield(TR, 'ConnectID')
            for jj = 1:numel(TR)
                TR(jj).ConnectID = TR(jj).TrajectoryID;
            end
        end
        FileList{dd} = [FileList{dd}(:); NewFileList{dd}(ii)];
        if (numel(TR) < 2)
            continue    % skip empty results (as in the original analysis)
        end
        DEstimator.TR = TR;
        DS = DEstimator.estimateDiffusionConstant();
        if ((DS(2).DiffusionConstant < DMinMax(1)) ...
                || (DS(2).DiffusionConstant > DMinMax(2)))
            continue    % excluded by D thresholds
        end
        KeptFiles{dd} = [KeptFiles{dd}; NewFileList{dd}(ii)];
        DAllCells{dd} = [DAllCells{dd}; ...
            [DS(2).DiffusionConstant, DS(2).DiffusionConstantSE]];
    end

    % Recompute condition-level results over all kept cells (old + new).
    TRAll = struct([]);
    for ii = 1:numel(KeptFiles{dd})
        load(KeptFiles{dd}{ii}, 'TR')
        if ~isfield(TR, 'ConnectID')
            for jj = 1:numel(TR)
                TR(jj).ConnectID = TR(jj).TrajectoryID;
            end
        end
        TRAll = smi_core.TrackingResults.catTR(TRAll, TR);
    end

    % MSD-based estimates.
    DEstimator.TR = TRAll;
    DS = DEstimator.estimateDiffusionConstant();
    MSDEnsemble{dd} = DEstimator.MSDEnsemble;
    DAllTrajNew = [DS(1).DiffusionConstant, DS(1).DiffusionConstantSE];
    KeepBool = ((DAllTrajNew(:, 1) >= DMinMax(1)) ...
        & (DAllTrajNew(:, 1) <= DMinMax(2)));
    DAllTraj{dd} = DAllTrajNew(KeepBool, :);
    NTraj(dd) = sum(KeepBool);
    MeanPhotonsAllTrajNew = ...
        arrayfun(@(i) mean(TRAll(i).Photons), 1:numel(TRAll));
    MeanPhotonsAllTraj{dd} = MeanPhotonsAllTrajNew(KeepBool);
    DEnsemble(dd) = DS(2).DiffusionConstant;
    DEnsembleSE(dd) = DS(2).DiffusionConstantSE;
    NCells(dd) = numel(KeptFiles{dd});
    NJumpsFit(dd) = sum(DEstimator.MSDEnsemble.NPoints(DEstimator.NFitPoints));

    % CDF-of-jumps estimate.
    DEstimatorCDF.TR = TRAll;
    DiffusionStruct{dd} = DEstimatorCDF.estimateDiffusionConstant();
    CDFEnsemble{dd} = DEstimatorCDF.MSDEnsemble;
    NTrajEnsemble(dd) = numel(TRAll);
    NJumps(dd) = numel(DEstimatorCDF.MSDEnsemble.CDFOfJumps);
end

% Save the updated data (same variables as before; DFit/DFitSE only if they
% exist in the loaded file -- re-run the frame-window section to refresh them).
SaveVars = {'FileList', 'KeptFiles', 'ConditionNames', 'NConditions', ...
    'BaseDir', 'experimentName', 'chName', 'plotTitle', 'SaveDir', ...
    'FrameRate', 'PixelSize', 'DMinMax', 'NComponents', 'SMF', ...
    'MSDEnsemble', ...
    'DAllTraj', 'DAllCells', 'MeanPhotonsAllTraj', ...
    'DEnsemble', 'DEnsembleSE', ...
    'NCells', 'NTraj', 'NTrajEnsemble', 'NJumps', 'NJumpsFit', ...
    'CDFEnsemble', 'DiffusionStruct'};
if exist('DFit', 'var') && exist('DFitSE', 'var')
    SaveVars = [SaveVars, {'DFit', 'DFitSE'}];
end
% Save as -v7; fall back to -v7.3 if a variable exceeds the 2 GB limit.
lastwarn('');
save(AnalysisFile, SaveVars{:});
if ~isempty(lastwarn)
    fprintf('-v7 save warned (%s); re-saving as -v7.3.\n', lastwarn);
    save(AnalysisFile, SaveVars{:}, '-v7.3');
end
% Also update the file list saved alongside the analysis.
save(fullfile(SaveDir, strcat(experimentName, "FileList")), "FileList");
fprintf('Saved updated diffusion analysis data to %s\n', AnalysisFile);
fprintf('Cells per condition now: %s\n', mat2str(NCells'));

%% Loop through the TR structures and estimate the diffusion constants for a
%% time (frame) window.

% calculate max frame number of all tractories
WindowWidth = 10;      % number of frames in sliding window
MaxTrajLength = 500;   % maximum number of frames to consider overall

SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FrameRate = FrameRate; % frames per second
SMF.Data.PixelSize = PixelSize; % um
DEst = smi_stat.DiffusionEstimator([], SMF);
DEst.UnitFlag = true;
DEst.FrameLagRange = [1, inf];
DEst.NFitPoints = 5;
DEst.FitMethod = 'WeightedLS';
DEst.FitTarget = 'MSD';
DEst.NComponents = 1;

% Estimated diffusion constant per time point and condition (DFit),
% and the SE produced in computing this value (DFitSE).
DFit   = zeros(MaxTrajLength, NConditions);
DFitSE = zeros(MaxTrajLength, NConditions);
for MinFrameNum = 1:MaxTrajLength
   fprintf(' [%d]', MinFrameNum);
   MaxFrameNum = MinFrameNum + WindowWidth - 1;
   for nn = 1:NConditions
      fprintf(' %d', nn);
      if nn == NConditions
         fprintf('\n');
      end
      % Load the TR structures for this condition.
      NTRCurrent = numel(FileList{nn});
      TRAll = struct([]);
      for ii = 1:NTRCurrent
         % Load the TR structure and modify some properties as appropriate.
         load(FileList{nn}{ii}, 'TR');
         if ~isfield(TR, 'ConnectID')
            for jj = 1:numel(TR)
               TR(jj).ConnectID = TR(jj).TrajectoryID;
            end   % jj
         end

         % Choose a frame number window from the TR structure.
         TR = smi_core.TrackingResults.windowTR(TR, MinFrameNum, MaxFrameNum);
         TR = smi_core.TrackingResults.windowStartTR(TR, MaxStartingFrameNum);

         % Store the current results in the concatenated arrays.
         %TRAll = smi_core.TrackingResults.catTR(TRAll, TR);
         for jj = 1:numel(TR)
            % Eliminate TRs that have zero jumps after windowing.
            if numel(TR(jj).FrameNum) > 0
               TRAll = smi_core.TrackingResults.catTR(TRAll, TR(jj));
            end
         end   % jj
      end   % ii [NTRCurrent]

      try
         % Find ensemble MSD
         [MSDSingleTrajtmp, MSDEnsembletmp] = ...
            smi_stat.DiffusionEstimator.computeMSD(TRAll);

         % do ensemble MSD fit to get D using 1-3 frame lags
         % MSDFit() = ....
         DEst.TR = TRAll;
         DiffusionStruct = DEst.estimateDiffusionConstant();
         % DiffusionStruct(1).Name = 'Trajectory'
         % DiffusionStruct(2).Name = 'Ensemble'
         %
         % Ensemble windowed diffusion constant (note DEst.FrameLagRange):
         DFit(MinFrameNum, nn)   = DiffusionStruct(2).DiffusionConstant;
         DFitSE(MinFrameNum, nn) = DiffusionStruct(2).DiffusionConstantSE;
      catch ME
         fprintf('\nDEst problem: MinFrameNum = %d, Condition = %d\n', ...
                 MinFrameNum, nn);
         fprintf('%s\n', ME.message);
         DFit(MinFrameNum, nn)   = NaN;
         DFitSE(MinFrameNum, nn) = NaN;
      end   % try
   end   % nn
end   % MinFrameNum

delta_t = 1 / SMF.Data.FrameRate;   % seconds per frame
DT = (1 : size(DFit, 1)) * delta_t;
PlotAxes = axes(figure());
hold(PlotAxes, 'on');
LegendEntries = cell(NConditions, 1);
ConditionColors = lines(NConditions);
for nn = 1 : NConditions
   DFit_nonan = DFit(:, nn);
   DFit_nonan = DFit_nonan(~isnan(DFit_nonan));
   fprintf('[%d] D = %f +- %f um^2/s\n', nn, ...
           mean(DFit_nonan), std(DFit_nonan));
   plot(DT', DFit(:, nn), '-', 'Color', ConditionColors(nn, :), ...
        'LineWidth', 2);
   LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
      ConditionNames{nn}, mean(DFit_nonan), std(DFit_nonan));
end
PlotAxes.FontWeight = 'bold';
legend(PlotAxes, LegendEntries, 'Location', 'best');
title( ...
   sprintf('Sliding Window (FrameLag = %d:%d, WindowWidth = %d)', ...
           DEst.FrameLagRange, WindowWidth));
xlabel('t (s)', 'Interpreter', 'tex', 'FontWeight', 'bold');
ylabel('D (\mu m^2 s^{-1})', 'Interpreter', 'tex', 'FontWeight', 'bold');
hold off
if ~isempty(SaveDir)
   saveas(gcf, fullfile(SaveDir, 'DSlidingWindow.png'));
   saveas(gcf, fullfile(SaveDir, 'DSlidingWindow.fig'));
   save(fullfile(SaveDir, 'DFit+SE.mat'), 'DFit', 'DFitSE');
end
