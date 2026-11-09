%% ExportCellDiffusionCSV
% Load a saved DiffusionAnalysisData.mat (from IMA_diffusion_script_A431_combined)
% and export the per-cell diffusion coefficients (the values plotted in
% DCoeffperCond) to CSV.
%
% Outputs (written next to the .mat file, or to OutDir if set):
%   <prefix>_CellDiffusion_long.csv  one row per cell:
%       Condition, ConditionIndex, CellIndex, D_um2_per_s, D_SE_um2_per_s, File
%   <prefix>_CellDiffusion_wide.csv  one column per condition (D only),
%       padded with empty cells; convenient for Prism/Excel.

%% Settings
AnalysisFile = '';   % full path to DiffusionAnalysisData.mat; '' = pick one
OutDir = '';         % '' = same folder as AnalysisFile

%% Load only the small variables we need (skips the large MSD structures).
if isempty(AnalysisFile)
    [fn, fp] = uigetfile('*.mat', 'Select DiffusionAnalysisData.mat');
    if isequal(fn, 0)
        error('No file selected.');
    end
    AnalysisFile = fullfile(fp, fn);
end
if isempty(OutDir)
    OutDir = fileparts(AnalysisFile);
end
S = load(AnalysisFile, 'DAllCells', 'KeptFiles', 'ConditionNames', ...
    'NConditions', 'experimentName', 'DAllTraj', 'MeanPhotonsAllTraj');
DAllCells = S.DAllCells;
DAllTraj = S.DAllTraj;
MeanPhotonsAllTraj = S.MeanPhotonsAllTraj;
KeptFiles = S.KeptFiles;
ConditionNames = S.ConditionNames;
NConditions = S.NConditions;
if isfield(S, 'experimentName') && ~isempty(S.experimentName)
    Prefix = char(S.experimentName);
else
    [~, Prefix] = fileparts(AnalysisFile);
end

%% Long format: one row per cell.
Condition = strings(0, 1);
ConditionIndex = zeros(0, 1);
CellIndex = zeros(0, 1);
D_um2_per_s = zeros(0, 1);
D_SE_um2_per_s = zeros(0, 1);
File = strings(0, 1);
for nn = 1:NConditions
    NCellsCurrent = size(DAllCells{nn}, 1);
    if NCellsCurrent ~= numel(KeptFiles{nn})
        warning(['Condition %d: %d D values but %d kept files; ', ...
            'file names may not line up.'], nn, NCellsCurrent, ...
            numel(KeptFiles{nn}));
    end
    Files = strings(NCellsCurrent, 1);
    NFiles = min(NCellsCurrent, numel(KeptFiles{nn}));
    Files(1:NFiles) = string(KeptFiles{nn}(1:NFiles));
    Condition = [Condition; repmat(string(ConditionNames{nn}), NCellsCurrent, 1)]; %#ok<AGROW>
    ConditionIndex = [ConditionIndex; repmat(nn, NCellsCurrent, 1)]; %#ok<AGROW>
    CellIndex = [CellIndex; (1:NCellsCurrent)']; %#ok<AGROW>
    D_um2_per_s = [D_um2_per_s; DAllCells{nn}(:, 1)]; %#ok<AGROW>
    D_SE_um2_per_s = [D_SE_um2_per_s; DAllCells{nn}(:, 2)]; %#ok<AGROW>
    File = [File; Files]; %#ok<AGROW>
end
LongTable = table(Condition, ConditionIndex, CellIndex, D_um2_per_s, ...
    D_SE_um2_per_s, File);
LongFile = fullfile(OutDir, sprintf('%s_CellDiffusion_long.csv', Prefix));
writetable(LongTable, LongFile);
fprintf('Wrote %d cells to %s\n', height(LongTable), LongFile);

%% Wide format: one column per condition (D only).
MaxCells = max(cellfun(@(X) size(X, 1), DAllCells));
WideData = nan(MaxCells, NConditions);
for nn = 1:NConditions
    WideData(1:size(DAllCells{nn}, 1), nn) = DAllCells{nn}(:, 1);
end
VarNames = matlab.lang.makeValidName(string(ConditionNames), ...
    'ReplacementStyle', 'delete');
VarNames = matlab.lang.makeUniqueStrings(VarNames);
WideTable = array2table(WideData, 'VariableNames', cellstr(VarNames));
WideFile = fullfile(OutDir, sprintf('%s_CellDiffusion_wide.csv', Prefix));
writetable(WideTable, WideFile);
fprintf('Wrote wide table to %s\n', WideFile);

%% Per-trajectory D (DAllTraj): long format, one row per trajectory.
% Trajectories are pooled across all cells in a condition (after the DMinMax
% filter), so there is no per-cell file column.
TCondition = strings(0, 1);
TConditionIndex = zeros(0, 1);
TrajIndex = zeros(0, 1);
TD = zeros(0, 1);
TDSE = zeros(0, 1);
TMeanPhotons = zeros(0, 1);
for nn = 1:NConditions
    NTrajCurrent = size(DAllTraj{nn}, 1);
    Photons = nan(NTrajCurrent, 1);
    if numel(MeanPhotonsAllTraj{nn}) == NTrajCurrent
        Photons = MeanPhotonsAllTraj{nn}(:);
    else
        warning('Condition %d: photon counts do not match trajectories.', nn);
    end
    TCondition = [TCondition; repmat(string(ConditionNames{nn}), NTrajCurrent, 1)]; %#ok<AGROW>
    TConditionIndex = [TConditionIndex; repmat(nn, NTrajCurrent, 1)]; %#ok<AGROW>
    TrajIndex = [TrajIndex; (1:NTrajCurrent)']; %#ok<AGROW>
    TD = [TD; DAllTraj{nn}(:, 1)]; %#ok<AGROW>
    TDSE = [TDSE; DAllTraj{nn}(:, 2)]; %#ok<AGROW>
    TMeanPhotons = [TMeanPhotons; Photons]; %#ok<AGROW>
end
TrajTable = table(TCondition, TConditionIndex, TrajIndex, TD, TDSE, ...
    TMeanPhotons, 'VariableNames', {'Condition', 'ConditionIndex', ...
    'TrajIndex', 'D_um2_per_s', 'D_SE_um2_per_s', 'MeanPhotons'});
TrajFile = fullfile(OutDir, sprintf('%s_TrajDiffusion_long.csv', Prefix));
writetable(TrajTable, TrajFile);
fprintf('Wrote %d trajectories to %s\n', height(TrajTable), TrajFile);
