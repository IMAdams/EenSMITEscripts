%% ReplotDiffusionFromSaved
% Load a saved DiffusionAnalysisData.mat (from IMA_diffusion_script_A431_combined)
% and reproduce the diffusion plots without re-running the analysis or
% touching the TR files.
%
% Plots reproduced (same names as the original script):
%   MSDvsDeltat            MSD curves, color by group / linestyle by sub-condition
%   MSDvsDeltatSHADED      MSD curves with SEM shading
%   CDFvsr2                CDF of squared jumps
%   DComponentFits         CDF with N-component Brownian fits
%   DCoeffperCond          per-cell D box plot + points
%   TrajDCoeffperCond      per-trajectory D spread
%   TrajDCoeffHist         histogram of -log10(D) per trajectory
%   TrajPhotonsperCond     per-trajectory mean photons
%   DSlidingWindow         sliding-window D (only if DFit was saved)
% Also rewrites summary.txt via IMA_exportD_est_Nvals.
%
% Requires smite on the path (smi_vis.plotCDFs, smi_stat.DiffusionEstimator)
% and plotSpread.

%% Settings
AnalysisFile = '';    % full path to DiffusionAnalysisData.mat; '' = pick one
OutDir = '';          % '' = <folder of AnalysisFile>/replots
SavePlots = true;     % write .png and .fig for each plot
NSubConditions = 3;   % conditions per group (as in the original script)
DYLimCells = [0, 0.16];  % y-limits of the per-cell D plot; [] = automatic
WindowWidth = 10;     % sliding-window width (frames) used to make DFit; only
                      % used in the DSlidingWindow title
SlidingFrameLagRange = [1, inf];  % FrameLagRange used for DFit (title only)
rng(1);               % reproducible jitter in the per-cell plot

%% Load.
if isempty(AnalysisFile)
    [fn, fp] = uigetfile('*.mat', 'Select DiffusionAnalysisData.mat');
    if isequal(fn, 0)
        error('No file selected.');
    end
    AnalysisFile = fullfile(fp, fn);
end
if isempty(OutDir)
    OutDir = fullfile(fileparts(AnalysisFile), 'replots');
end
if SavePlots && ~isfolder(OutDir)
    mkdir(OutDir);
end
% Note: this deliberately ignores the SaveDir stored in the file, which may
% point at a drive/path that no longer exists.
load(AnalysisFile);
fprintf('Loaded %s\n', AnalysisFile);

PixelSizeUm = PixelSize;       % um/pixel
FrameRateHz = FrameRate;       % frames/s
NGroups = ceil(NConditions / NSubConditions);
ConditionColors = lines(NConditions);

DAllCellsIsolated = cellfun(@(X) X(:, 1), DAllCells, 'UniformOutput', false);
DAllTrajIsolated = cellfun(@(X) X(:, 1), DAllTraj, 'UniformOutput', false);

%% MSD vs. delta t (group colors, sub-condition linestyles).
Fig = figure();
PlotAxes = axes(Fig);
hold(PlotAxes, 'on')
LegendEntries = cell(NConditions, 1);
GroupColors = lines(NGroups);
LineStyles = {'-', '--', ':'};
for nn = 1:NConditions
    GroupIndex = ceil(nn / NSubConditions);
    SubIndex = mod(nn - 1, NSubConditions) + 1;
    plot(PlotAxes, MSDEnsemble{nn}.FrameLags / FrameRateHz, ...
        MSDEnsemble{nn}.MSD * (PixelSizeUm^2), ...
        'LineWidth', 2, ...
        'Color', GroupColors(GroupIndex, :), ...
        'LineStyle', LineStyles{SubIndex})
    LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
        ConditionNames{nn}, DEnsemble(nn), 1.96 * DEnsembleSE(nn));
end
title(PlotAxes, strcat(plotTitle, ' MSD'))
PlotAxes.FontWeight = 'bold';
PlotAxes.XLim = [0, 1];
legend(PlotAxes, LegendEntries, 'Location', 'northwest')
xlabel(PlotAxes, '\Deltat (s)', 'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxes, 'MSD (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
savePlot(Fig, OutDir, 'MSDvsDeltat', SavePlots);

%% MSD vs. delta t with SEM shading.
HasRawJumps = isfield(MSDEnsemble{1}, 'SquaredDisplacement') ...
    && isfield(MSDEnsemble{1}, 'FrameLagsAll');
if HasRawJumps
    Fig = figure();
    PlotAxes = axes(Fig);
    hold(PlotAxes, 'on')
    LegendEntries = cell(NConditions, 1);
    for nn = 1:NConditions
        tLags = MSDEnsemble{nn}.FrameLags / FrameRateHz;
        MSDPlot = MSDEnsemble{nn}.MSD * (PixelSizeUm^2);

        % SEM at each lag from the raw squared displacements.
        NLags = numel(MSDEnsemble{nn}.FrameLags);
        SEM = zeros(NLags, 1);
        for ii = 1:NLags
            SD = MSDEnsemble{nn}.SquaredDisplacement( ...
                MSDEnsemble{nn}.FrameLagsAll == MSDEnsemble{nn}.FrameLags(ii));
            SEM(ii) = std(SD) / sqrt(numel(SD));
        end
        SEMPlot = SEM * (PixelSizeUm^2);

        XPatch = [tLags(:); flip(tLags(:))];
        YPatch = [MSDPlot(:) + SEMPlot; flip(MSDPlot(:) - SEMPlot)];
        patch(PlotAxes, XPatch, YPatch, ConditionColors(nn, :), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(PlotAxes, tLags, MSDPlot, ...
            'LineWidth', 2, 'Color', ConditionColors(nn, :))
        LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
            ConditionNames{nn}, DEnsemble(nn), 1.96 * DEnsembleSE(nn));
    end
    PlotAxes.FontWeight = 'bold';
    PlotAxes.XLim = [0, 1];
    legend(PlotAxes, LegendEntries, 'Location', 'northwest')
    xlabel(PlotAxes, '\Deltat (s)', 'Interpreter', 'tex', 'FontWeight', 'bold')
    ylabel(PlotAxes, 'MSD (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
    savePlot(Fig, OutDir, 'MSDvsDeltatSHADED', SavePlots);
else
    warning(['MSDEnsemble has no SquaredDisplacement/FrameLagsAll; ', ...
        'skipping the SEM-shaded MSD plot.']);
end

%% CDFs of squared jumps.
SqDisp = cell(NConditions, 1);
CDFs = cell(NConditions, 1);
for nn = 1:NConditions
    SqDisp{nn} = CDFEnsemble{nn}.SortedSquaredDisp * (PixelSizeUm^2);
    if isfield(CDFEnsemble{nn}, 'CDFOfJumps')
        CDFs{nn} = CDFEnsemble{nn}.CDFOfJumps;
    else
        % CDFOfJumps is just the empirical CDF of the sorted jumps.
        NJ = numel(CDFEnsemble{nn}.SortedSquaredDisp);
        CDFs{nn} = (1:NJ)' / NJ;
    end
end
Fig = figure();
PlotAxes = axes(Fig);
hold(PlotAxes, 'on')
smi_vis.plotCDFs(PlotAxes, SqDisp, CDFs, ConditionColors)
PlotAxes.FontWeight = 'bold';
axis(PlotAxes, 'tight')
PlotAxes.XScale = 'log';
legend(PlotAxes, ConditionNames, 'Location', 'best')
xlabel(PlotAxes, 'r^2 (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxes, 'CDF(r^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
title(PlotAxes, strcat(plotTitle, ' CDF'))
savePlot(Fig, OutDir, 'CDFvsr2', SavePlots);

%% CDFs with component fits.
Fig = figure();
PlotAxes = axes(Fig);
hold(PlotAxes, 'on')
ModelCDFs = cell(NConditions, 1);
LegendEntries = cell(2 * NConditions, 1);
LegendEntries(1:NConditions) = ConditionNames;
for nn = 1:NConditions
    ModelCDFs{nn} = smi_stat.DiffusionEstimator.brownianJumpCDF(...
        [DiffusionStruct{nn}(2).DiffusionConstant, ...
        DiffusionStruct{nn}(2).PopulationRatios], ...
        CDFEnsemble{nn}.SortedSquaredDisp * (PixelSizeUm^2), ...
        CDFEnsemble{nn}.FrameLags / FrameRateHz, ...
        CDFEnsemble{nn}.NPoints, ...
        CDFEnsemble{nn}.LocVarianceSum * (PixelSizeUm^2));
    LegendEntries{nn + NConditions} = ...
        [sprintf('%s fit', ConditionNames{nn}), newline, ...
        sprintf('D_1 = %.4g, D_2 = %.4g, D_3 = %.4g, Pop_1 = %.4f, Pop_2 = %.4f, Pop_3 = %.4f', ...
        DiffusionStruct{nn}(2).DiffusionConstant, ...
        DiffusionStruct{nn}(2).PopulationRatios)];
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
savePlot(Fig, OutDir, 'DComponentFits', SavePlots);

%% Per-cell diffusion coefficients.
Fig = figure();
PlotAxes = axes(Fig);
hold(PlotAxes, 'on')
AllD = vertcat(DAllCellsIsolated{:});
GroupLabels = [];
for nn = 1:NConditions
    GroupLabels = [GroupLabels; repmat(nn, numel(DAllCellsIsolated{nn}), 1)]; %#ok<AGROW>
end
boxplot(PlotAxes, AllD, GroupLabels, ...
    'Labels', ConditionNames, ...
    'Colors', 'k', ...
    'Symbol', '');
for nn = 1:NConditions
    xJitter = 0.15 * (rand(size(DAllCellsIsolated{nn})) - 0.5);
    scatter(PlotAxes, nn + xJitter, DAllCellsIsolated{nn}, ...
        45, ConditionColors(ceil(nn / NSubConditions), :), 'filled', ...
        'MarkerFaceAlpha', 0.7)
end
PlotAxes.XLim = [0, NConditions + 1];
PlotAxes.XTick = 1:NConditions;
PlotAxes.XTickLabel = ConditionNames;
PlotAxes.XTickLabelRotation = 45;
PlotAxes.FontWeight = 'bold';
if ~isempty(DYLimCells)
    PlotAxes.YLim = DYLimCells;
end
PlotAxes.FontSize = 15;
ylabel(PlotAxes, 'D (\mum^2s^{-1})', 'Interpreter', 'tex', ...
    'FontWeight', 'bold')
title(PlotAxes, strcat(plotTitle, ' Cell diffusion coefficients from MSD fits'))
savePlot(Fig, OutDir, 'DCoeffperCond', SavePlots);

%% Per-trajectory diffusion coefficients.
Fig = figure();
PlotAxes = axes(Fig);
hold(PlotAxes, 'on')
plotSpread(PlotAxes, DAllTrajIsolated, ...
    'spreadWidth', 0.5, 'distributionColors', ConditionColors)
for nn = 1:NConditions
    plot(PlotAxes, nn + 0.25 * [-1, 1], mean(DAllTrajIsolated{nn}) * [1, 1], ...
        'Color', [0, 0, 0], 'LineWidth', 2)
end
PlotAxes.XLim = [0, NConditions + 1];
PlotAxes.XTick = 1:NConditions;
PlotAxes.XTickLabel = ConditionNames;
PlotAxes.XTickLabelRotation = 45;
PlotAxes.FontWeight = 'bold';
ylabel(PlotAxes, 'D (\mum^2s^{-1})', 'Interpreter', 'tex', 'FontWeight', 'bold')
title(PlotAxes, 'Trajectory diffusion coefficients from MSD fits')
savePlot(Fig, OutDir, 'TrajDCoeffperCond', SavePlots);

%% Histogram of -log10(D) per trajectory.
Fig = figure();
PlotAxesHist = axes(Fig);
hold(PlotAxesHist, 'on')
DLogAll = cellfun(@(X) -log10(X), DAllTrajIsolated, 'UniformOutput', false);
AllDLog = vertcat(DLogAll{:});
AllDLog = AllDLog(isfinite(AllDLog));   % D = 0 gives Inf
BinEdges = linspace(min(AllDLog), max(AllDLog), 30);
for nn = 1:NConditions
    histogram(PlotAxesHist, DLogAll{nn}, BinEdges, ...
        'FaceColor', ConditionColors(nn, :), ...
        'FaceAlpha', 0.5, ...
        'EdgeColor', 'none', ...
        'Normalization', 'probability', ...
        'DisplayName', ConditionNames{nn})
end
xlabel(PlotAxesHist, '-log_{10}(D) (-log_{10}(\mum^2s^{-1}))', ...
    'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxesHist, 'probability', 'FontWeight', 'bold')
title(PlotAxesHist, 'D from MSD fits ')
legend(PlotAxesHist, 'Location', 'best')
PlotAxesHist.FontWeight = 'bold';
% Secondary x-axis showing actual D values.
ax2 = axes(Fig, 'Position', PlotAxesHist.Position, ...
    'XAxisLocation', 'top', ...
    'Color', 'none', ...
    'YTick', []);
ax2.XLim = PlotAxesHist.XLim;
TickVals = PlotAxesHist.XTick;
ax2.XTick = TickVals;
ax2.XTickLabel = arrayfun(@(x) sprintf('%.3g', 10^(-x)), TickVals, ...
    'UniformOutput', false);
ax2.FontWeight = 'bold';
xlabel(ax2, 'D (\mum^2s^{-1})', 'Interpreter', 'tex', 'FontWeight', 'bold')
savePlot(Fig, OutDir, 'TrajDCoeffHist', SavePlots);

%% Per-trajectory mean photons.
Fig = figure();
PlotAxes = axes(Fig);
hold(PlotAxes, 'on')
plotSpread(PlotAxes, MeanPhotonsAllTraj, ...
    'spreadWidth', 0.5, 'distributionColors', ConditionColors)
for nn = 1:NConditions
    plot(PlotAxes, nn + 0.25 * [-1, 1], mean(MeanPhotonsAllTraj{nn}) * [1, 1], ...
        'Color', [0, 0, 0], 'LineWidth', 2)
end
PlotAxes.XLim = [0, NConditions + 1];
PlotAxes.XTick = 1:NConditions;
PlotAxes.XTickLabel = ConditionNames;
PlotAxes.XTickLabelRotation = 45;
PlotAxes.FontWeight = 'bold';
ylabel(PlotAxes, 'mean photons', 'Interpreter', 'tex', 'FontWeight', 'bold')
title(PlotAxes, 'Trajectory mean photons')
savePlot(Fig, OutDir, 'TrajPhotonsperCond', SavePlots);

%% Sliding-window D (only if it was computed and saved).
if exist('DFit', 'var') && exist('DFitSE', 'var')
    DT = (1:size(DFit, 1)) / FrameRateHz;   % seconds
    Fig = figure();
    PlotAxes = axes(Fig);
    hold(PlotAxes, 'on');
    LegendEntries = cell(NConditions, 1);
    for nn = 1:NConditions
        DFitNoNaN = DFit(~isnan(DFit(:, nn)), nn);
        plot(PlotAxes, DT', DFit(:, nn), '-', ...
            'Color', ConditionColors(nn, :), 'LineWidth', 2);
        LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
            ConditionNames{nn}, mean(DFitNoNaN), std(DFitNoNaN));
    end
    PlotAxes.FontWeight = 'bold';
    legend(PlotAxes, LegendEntries, 'Location', 'best');
    title(PlotAxes, sprintf( ...
        'Sliding Window (FrameLag = %d:%d, WindowWidth = %d)', ...
        SlidingFrameLagRange, WindowWidth));
    xlabel(PlotAxes, 't (s)', 'Interpreter', 'tex', 'FontWeight', 'bold');
    ylabel(PlotAxes, 'D (\mu m^2 s^{-1})', 'Interpreter', 'tex', ...
        'FontWeight', 'bold');
    savePlot(Fig, OutDir, 'DSlidingWindow', SavePlots);
else
    fprintf('No DFit/DFitSE in the file; skipping the sliding-window plot.\n');
end

%% Summary of N values (same helper as the original script).
if SavePlots
    try
        IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, ...
            DAllTrajIsolated, fullfile(OutDir, 'summary.txt'), ...
            'Format', 'detailed')
    catch ME
        warning('Could not write summary.txt: %s', ME.message);
    end
end
fprintf('Done. Plots written to %s\n', OutDir);

%% Local functions
function savePlot(Fig, OutDir, Name, SavePlots)
    if SavePlots
        saveas(Fig, fullfile(OutDir, [Name, '.png']));
        saveas(Fig, fullfile(OutDir, [Name, '.fig']));
    end
end
