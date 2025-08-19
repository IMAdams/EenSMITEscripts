 % This script is used to estimate diffusion constants from TR structures.

%% Define some directories and parameters.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FrameRate = 20;
SMF.Data.PixelSize = 0.1667;
DMinMax = [0, inf]; % um^2/s, [min., max.] for MSD results (cell and trajectory)
NComponents = 3; % number of diffusing components for CDF/MLE fitting
BaseDir = 'C:\Users\imadams\Documents\smite workspace\20250702_CHO-HA-EGFR_CV_dEGF\Results0729';
PlotSaveDir = 'C:\Users\imadams\Documents\smite workspace\20250702_CHO-HA-EGFR_CV_dEGF\Results0729\D';
% 
expName = '0702_CHOHAEGFR';

ConditionNames = {...
        'HA-EGFR-WT Resting'...
        'HA-EGFR-WT +EGF'...
        'HA-EGFR-L858R Resting'...
        'HA-EGFR-L858R +EGF'};
    % 'Cell 1'...
    % 'Cell 2'...
    % 'Cell 3'...
    % 'Cell 4'...
    % 'Cell 5'...
    % 'Cell 7'...
    % 'Cell 8'...
    % 'Cell 9'...
    % 'Cell 10'...
    % 'Cell 11'...
    % 'Cell 12'...
    % 'Cell 13'...
    % 'Cell 14'...
    % 'Cell 15'...
    % 'Cell 16'...
    % 'Cell 17'...
    % 'Cell 18'...
    % 'Cell 19'...
    % 'Cell 20'...
    % {'Anti-ALFA Resting'...
    % 'Anti-ALFA +EGF'...
    % 'Anti-ALFA Resting'...
    % 'Anti-ALFA +TGF-B'};
    % {'CHO-HA-EGFR Resting'...
    % 'CHO-HA-EGFR +EGF'...
    % 'CHO-HA-EGFR-L858R Resting'...
    % 'CHO-HA-EGFR-L858R +EGF'};
   % {'HA-EzrinABD'...
    %    'HA-EzrinABD+LatB'}
    % {'HER2-WT Resting'...
  %  {'HER2-WT +EGF'...
    % 'HER2-S310F Resting'...
  %   'HER2-S310F +EGF'};
    % {'RBL-HA-EzrinABD'...
    % 'RBL-HA-EzrinABD + LatB'};
% {'HER2-WT EGF-QD'...
% 'HER2-S310F EGF-QD'};
% ConditionNames = ...
%     {'CHOK1 HA-EGFR-WT'...
%     'CHOK1 HA-EGFR-WT + EGF'...
%     'CHOK1 HA-EGFR-L858R'...
%     'CHOK1 HA-EGFR-L858R + EGF'};

filterBy = '1_Results.mat';

NConditions = numel(ConditionNames);
FileList = cell(NConditions, 1);
for ii = 1:NConditions
    fprintf('Pick files for condition %s\n', ConditionNames{ii})
    FileList{ii} = uipickfiles('FilterSpec', BaseDir, ...
        'REFilter', filterBy, ...
        'Prompt', ...
        sprintf('Pick files for condition %s\n', ConditionNames{ii}));
end

Verbose = 1;

%% Loop through the TR structures and estimate the diffusion constants.
DEstimator = smi_stat.DiffusionEstimator([], SMF, 3);
DEstimator.UnitFlag = true;
DEstimator.FrameLagRange = [1, 50];
DEstimator.NFitPoints = 5;
DEstimator.FitMethod = 'WeightedLS'; 
DEstimator.FitTarget =  'MSD';
DEstimator.Verbose = Verbose;
DEstimator.EstimateSEs = 1;
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
        % Cell_id = IMA_Cell_id(string(FileList{dd}(ii)));
        % if ~isfield(TR, 'cell_id')
        %     for jj = 1:numel(TR)
        %         TR(jj).cell_id = Cell_id;
        %     end
        % end
        % Estimate the trajectory-wise diffusion constants (one diffusion
        % constant per trajectory, excluding trajectories that were too

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
    DEstimator.NComponents = 1; % must be 1 for MSD fitting
    % DiffusionStruct{dd} = DEstimator.estimateDiffusionConstant();
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
% 
% Repeat above but estimate D from CDF of jumps (instead of fitting an MSD)
DEstimator = smi_stat.DiffusionEstimator([], SMF);
DEstimator.UnitFlag = true;
DEstimator.FrameLagRange = [5,5];
DEstimator.NFitPoints = inf;
% DEstimator.FitMethod = 'LS';
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
        % Create a field for the cell ID
        Cell_id = IMA_Cell_id(string(FileList{dd}(ii)));
        if ~isfield(TR, 'ConnectID')
            for jj = 1:numel(TR)
                TR(jj).ConnectID = TR(jj).TrajectoryID;
            end
        end
        % if ~isfield(TR, 'cell_id')
        %     for jj = 1:numel(TR)
        %         TR(jj).cell_id = Cell_id;
        %     end
        % end
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


% Plot some results.
% Plot the MSD for all conditions (ensemble over trajectories and cells).
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
LegendEntries = cell(NConditions, 1);
ConditionColors = lines(NConditions);
for nn = 1:NConditions
    plot(PlotAxes, MSDEnsemble{nn}.FrameLags/SMF.Data.FrameRate, ...
        MSDEnsemble{nn}.MSD*(SMF.Data.PixelSize^2), ...
        'LineWidth', 2, 'Color', ConditionColors(nn, :))
    LegendEntries{nn} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
        ConditionNames{nn}, DEnsemble(nn), 1.96*DEnsembleSE(nn));
end
PlotAxes.FontWeight = 'bold';
PlotAxes.XLim = [0, 2];
legend(PlotAxes, LegendEntries, 'Location', 'best')
xlabel(PlotAxes, '\Deltat (s)', 'Interpreter', 'tex', 'FontWeight', 'bold')
ylabel(PlotAxes, 'MSD (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
% Save the plot as .png and .fig
saveas(gcf, fullfile(PlotSaveDir, 'MeanSquaredDisplacement.png'))
savefig(gcf, fullfile(PlotSaveDir,'MeanSquaredDisplacement.fig'))

% Plot the CDFs for all conditions (ensemble over trajectories and cells).
SqDisp = cell(NConditions, 1);
CDFs = cell(NConditions, 1);
for nn = 1:NConditions
    SqDisp{nn} = CDFEnsemble{nn}.SortedSquaredDisp * (SMF.Data.PixelSize^2);
    CDFs{nn} = CDFEnsemble{nn}.CDFOfJumps;
end

% export for biostats

% IMA_extractCDFJumpsToh5(CDFEnsemble, char(fullfile(PlotSaveDir, strcat(expName, string(nn), ".h5"))))

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
% Save the plot as .png and .fig
saveas(gcf, fullfile(PlotSaveDir, 'CDF.png'))
savefig(gcf, fullfile(PlotSaveDir,'CDF.fig'))

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
    % LegendEntries{nn + NConditions} = sprintf('%s fit', ConditionNames{nn});
end
writecell(LegendEntries, fullfile(PlotSaveDir, "D_legend.txt"));

% IMA_exportStructToText(LegendEntries, fullfile(PlotSaveDir, "D_legend.txt"));
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
% Save the plot as .png and .fig
saveas(gcf, fullfile(PlotSaveDir, 'CDFfit.png'))
savefig(gcf, fullfile(PlotSaveDir,'CDFfit.fig'))

% Plot the diffusion constants for each cell.
PlotAxes = axes(figure());
hold(PlotAxes, 'on')
DAllCellsIsolated = cellfun(@(X) X(:, 1), DAllCells, ...
    'UniformOutput', false); % isolate D only, not D_SE
plotSpread(PlotAxes, DAllCellsIsolated, ...
    'spreadWidth', 0.5, 'distributionColors', ConditionColors)
for nn = 1:NConditions
    % Add a horizontal line at the mean value.
    plot(PlotAxes, nn + 0.25*[-1, 1], mean(DAllCellsIsolated{nn})*[1, 1], ...
        'Color', [0, 0, 0], 'LineWidth', 2)
end
PlotAxes.XLim = [0, NConditions+1];
PlotAxes.XTick = 1:NConditions;
PlotAxes.XTickLabel = ConditionNames;
PlotAxes.XTickLabelRotation = 45;
PlotAxes.FontWeight = 'bold';
ylabel(PlotAxes, 'D (\mum^2s^{-1})', 'Interpreter', 'tex', ...
    'FontWeight', 'bold')
title(PlotAxes, 'Cell diffusion coefficients from MSD fits')
% Save the plot as .png and .fig
saveas(gcf, fullfile(PlotSaveDir, 'CellD.png'))
savefig(gcf, fullfile(PlotSaveDir,'CellD.fig'))

% Plot the diffusion constants for each trajectory.
%
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
saveas(gcf, fullfile(PlotSaveDir, 'D_MSD.png'))
savefig(gcf, fullfile(PlotSaveDir,'D_MSD.fig'))

%
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
saveas(gcf, fullfile(PlotSaveDir, 'photonsTraj.png'))
savefig(gcf, fullfile(PlotSaveDir,'photonsTraj'))
for nn = 1:NConditions
   figure();
   hold on
   plot(MeanPhotonsAllTraj{nn}, DAllTrajIsolated{nn}, 'k.');
   title(ConditionNames{nn});
   xlabel('mean photons per trajectory');
   ylabel('D (\mum^2s^{-1}) per trajectory');
   hold off
end
% save values

outputFilename = strcat('N_numbers',expName,'.txt');
outputFilename = fullfile(PlotSaveDir, outputFilename );
IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, outputFilename)
%% test MLE
% Prepare the DiffusionEstimator class and estimate D.
% NOTE: For LikelihoodOfJumps maximization, the important parameters are
%       'DiffusionModel' and 'FrameLagRange'.
% SMF = smi_core.SingleMoleculeFitting; % used for PixelSize and FrameRate
DE = smi_stat.DiffusionEstimator(TR, SMF);
DE.FitTarget = 'LikelihoodOfJumps';
DE.DiffusionModel = 'Brownian';
DE.NComponents = 2;
DE.FrameLagRange = [1, 5]; % range of frame lags computed in MSD
DE.EstimateSEs = true; % estimate standard errors of fit D values
DE.FitIndividualTrajectories = false; % fit each trajectory MSD
DE.UnitFlag = false;
DE.Verbose = 1;
DE.estimateDiffusionConstant();

% 
% PlotAxes = axes(figure());
% hold(PlotAxes, 'on')
% LegendEntries = cell(NConditions, 1);
% ConditionColors = lines(NConditions);
% plot(PlotAxes, DE.MSDEnsemble.FrameLags/SMF.Data.FrameRate, ...
%         DE.MSDEnsemble.MSD*(SMF.Data.PixelSize^2), ...
%         'LineWidth', 2, 'Color', ConditionColors(1, :))
%     LegendEntries{1} = sprintf('%s (D=%.4f\\pm%.4f \\mum^2s^{-1})', ...
%         ConditionNames{1}, DEnsemble(1), 1.96*DEnsembleSE(1));
% end
% PlotAxes.FontWeight = 'bold';
% PlotAxes.XLim = [0, 5];
% legend(PlotAxes, LegendEntries, 'Location', 'best')
% xlabel(PlotAxes, '\Deltat (s)', 'Interpreter', 'tex', 'FontWeight', 'bold')
% ylabel(PlotAxes, 'MSD (\mum^2)', 'Interpreter', 'tex', 'FontWeight', 'bold')
% % Save the plot as .png and .fig
% saveas(gcf, fullfile(PlotSaveDir, 'MeanSquaredDisplacement.png'))
% savefig(gcf, fullfile(PlotSaveDir,'MeanSquaredDisplacement.fig'))