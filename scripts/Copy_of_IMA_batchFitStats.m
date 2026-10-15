%% Trajectory stats across experimental conditions
% Picks files per condition, aggregates trajectory stats from SMD structure(), tracks
% cell/file provenance, plots per-condition QC figures, prints summary
% stats, and produces a cross-condition comparison figure.

%% --- Config ---------------------------------------------------------
BaseDir = 'O:\Cell Path\Lidke Lab\IMAdams\Data\HeLa-ALFA-EGFR-KI\SPT\20260729_HeLaALFAEGFR_HATBR2\01_data_dot_mat\Results_SMF070726';
experimentName = 'ch2_EGF_TGFB1';
chName = '2'; % Qdot-655 -> '2'; Qdot-585 -> '1'
ConditionNames = ...
   {'ALFA-655 Rest' ...
    'ALFA-655 +EGF' ...
    'ALFA-655 +TGF\beta1' ...
   };
NConditions = numel(ConditionNames);
SaveDir = fullfile(BaseDir, strcat('D-analysis_', experimentName, '_ch', chName));
filterBy = strcat(chName, '_Results.mat');

% Fields to pull from each SMD and their display titles.
FieldNames = {'PValue', 'Photons', 'Bg', 'X_SE', 'Y_SE'};
FieldTitles = {'P-value', 'Photons', 'Background', 'X localization SE', 'Y localization SE'};

%% --- Run the file picker per condition -------------------------------
FileList = cell(NConditions, 1);
for ii = 1:NConditions
    fprintf('Pick files for condition %s\n', ConditionNames{ii})
    FileList{ii} = uipickfiles('FilterSpec', BaseDir, ...
        'REFilter', filterBy, ...
        'Prompt', ...
        sprintf('Pick files for condition %s\n', ConditionNames{ii}));
end

% If SaveDir is non-empty but doesn't exist yet, create the directory.
if ~isempty(SaveDir) && ~isfolder(SaveDir)
    mkdir(SaveDir);
end

% Keep a record of exactly which files went into which condition.
save(fullfile(SaveDir, [experimentName '_FileList.mat']), 'FileList', 'ConditionNames', 'BaseDir', 'filterBy');

%% --- Aggregate fit stats per condition, tracking file/cell ID --------
CondStats = struct();
for ii = 1:NConditions
    files = FileList{ii};
    nFiles = numel(files);

    aggregated = struct();
    for f = 1:numel(FieldNames)
        aggregated.(FieldNames{f}) = [];
    end
    aggregated.CellID = [];

    for k = 1:nFiles
        S = load(files{k});
        n = numel(S.SMD.(FieldNames{1}));
        for f = 1:numel(FieldNames)
            aggregated.(FieldNames{f}) = [aggregated.(FieldNames{f}); S.SMD.(FieldNames{f})(:)];
        end
        aggregated.CellID = [aggregated.CellID; repmat(k, n, 1)];
    end

    CondStats.(matlab.lang.makeValidName(ConditionNames{ii})) = aggregated;
end

%% --- Per-condition QC figure (subplots, log-scale P-value) -----------
for ii = 1:NConditions
    condField = matlab.lang.makeValidName(ConditionNames{ii});
    data = CondStats.(condField);

    fig = figure('Name', ['FitStats_' ConditionNames{ii}]);
    for f = 1:numel(FieldNames)
        subplot(2, 3, f);
        histogram(data.(FieldNames{f}));
        title(FieldTitles{f});
        if strcmp(FieldNames{f}, 'PValue')
            set(gca, 'YScale', 'log');
        end
    end
    sgtitle(sprintf('%s (n = %d localizations, %d files)', ...
        ConditionNames{ii}, numel(data.PValue), numel(FileList{ii})), 'Interpreter', 'tex');

    figName = sprintf('%s_ch%s_FitStats_%s', experimentName, chName, condField);
    savefig(fig, fullfile(SaveDir, [figName '.fig']));
    saveas(fig, fullfile(SaveDir, [figName '.png']));
end

%% --- Summary stats (median / IQR) per condition ----------------------
fprintf('\n=== Fit stat summary (median [IQR]) ===\n');
for ii = 1:NConditions
    condField = matlab.lang.makeValidName(ConditionNames{ii});
    data = CondStats.(condField);
    fprintf('\n%s:\n', ConditionNames{ii});
    for f = 1:numel(FieldNames)
        vals = data.(FieldNames{f});
        med = median(vals, 'omitnan');
        iqrVal = iqr(vals);
        fprintf('  %-20s median = %.4g   IQR = %.4g\n', FieldTitles{f}, med, iqrVal);
    end
end

%% --- Cross-condition comparison figure (boxplots + violin for P-value) --
fig = figure('Name', [experimentName '_ConditionComparison']);
for f = 1:numel(FieldNames)
    subplot(2, 3, f);

    allVals = [];
    allGroups = [];
    for ii = 1:NConditions
        condField = matlab.lang.makeValidName(ConditionNames{ii});
        vals = CondStats.(condField).(FieldNames{f});
        allVals = [allVals; vals];
        allGroups = [allGroups; repmat(ii, numel(vals), 1)];
    end

    if strcmp(FieldNames{f}, 'PValue')
        % P-value is bounded on [0,1] and often piles up near 0 -- a
        % violin shows that shape far better than a boxplot.
        localViolin(allVals, allGroups, ConditionNames, [0 1]);
        set(gca, 'YScale', 'log');
    else
        % 'Symbol','' suppresses the individual outlier markers -- with
        % this many points per box they just clutter the plot; the
        % whiskers/IQR already convey the spread.
        boxplot(allVals, allGroups, 'Labels', ConditionNames, 'Symbol', '');
    end
    title(FieldTitles{f});
    xtickangle(20);
end
sgtitle(sprintf('%s: fit stats by condition', experimentName), 'Interpreter', 'none');

compFigName = sprintf('%s_ch%s_ConditionComparison', experimentName, chName);
savefig(fig, fullfile(SaveDir, [compFigName '.fig']));
saveas(fig, fullfile(SaveDir, [compFigName '.png']));

%% --- Save aggregated stats for later reuse ---------------------------
save(fullfile(SaveDir, [experimentName '_ch' chName '_FitStatsAggregated.mat']), ...
    'CondStats', 'ConditionNames', 'FieldNames', 'FieldTitles');

fprintf('\nDone. Figures and aggregated stats saved to:\n  %s\n', SaveDir);

%% --- Local functions --------------------------------------------------
function localViolin(vals, groups, groupNames, support)
%localViolin Draw a simple violin plot (kernel density, mirrored) per
%group on the current axes. vals/groups are paired column vectors;
%groupNames is a cell array of labels; support is an optional [lo hi]
%bound passed to ksdensity (e.g. [0 1] for a P-value, which is bounded
%and often piles up near 0 -- ksdensity's boundary correction keeps
%that pileup from leaking density below 0).
    if nargin < 4
        support = 'unbounded';
    elseif isnumeric(support)
        % Pad slightly so the support strictly contains the data.
        pad = 1e-3 * range(support);
        support = [support(1) - pad, support(2) + pad];
    end

    nGroups = numel(groupNames);
    maxHalfWidth = 0.35;
    hold on;
    for g = 1:nGroups
        gVals = vals(groups == g);
        gVals = gVals(~isnan(gVals));
        if numel(gVals) < 2
            continue;
        end

        if isnumeric(support)
            [f, xi] = ksdensity(gVals, 'Support', support, ...
                'BoundaryCorrection', 'reflection');
        else
            [f, xi] = ksdensity(gVals);
        end
        f = f / max(f) * maxHalfWidth;

        patch([g - f, fliplr(g + f)], [xi, fliplr(xi)], [0.30 0.60 0.90], ...
            'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 0.75);

        medVal = median(gVals);
        plot([g - maxHalfWidth, g + maxHalfWidth], [medVal medVal], ...
            'r-', 'LineWidth', 1.5);
    end
    hold off;

    set(gca, 'XTick', 1:nGroups, 'XTickLabel', groupNames);
    xlim([0.5, nGroups + 0.5]);
end