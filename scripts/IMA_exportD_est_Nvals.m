function IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, MSDfitTraj, outputFilename, varargin)
% EXPORTCDFSUMMARYSTATS Export CDF and MSD-fit summary statistics to text file
%
% Usage:
%   IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, MSDfitTraj, 'summary.txt')
%   IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, MSDfitTraj, 'summary.txt', 'Format', 'table')
%   IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, MSDfitTraj, 'summary.txt', 'Format', 'detailed')
%
% Inputs:
%   ConditionNames  - 1×n cell array of condition names
%   NCells          - n×1 vector of number of cells per condition
%   NTraj           - n×1 vector of number of trajectories per condition
%   CDFs            - n×1 cell array of CDF data structures
%   MSDfitTraj      - n×1 cell array of trajectory-wise diffusion coefficients
%                     derived from MSD fit (each element is a numeric vector
%                     of D estimates per trajectory for that condition)
%   outputFilename  - String, name of output text file
%
% Optional Parameters:
%   'Format'        - 'table' (default): tabular format
%                   - 'detailed': detailed format with descriptions
%                   - 'compact': minimal format
%   'Delimiter'     - Delimiter for table format (default: '\t' for tab)

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'ConditionNames', @iscell);
    addRequired(p, 'NCells', @isnumeric);
    addRequired(p, 'NTraj', @isnumeric);
    addRequired(p, 'CDFs', @iscell);
    addRequired(p, 'MSDfitTraj', @iscell);
    addRequired(p, 'outputFilename', @ischar);
    addParameter(p, 'Format', 'table', @(x) ismember(x, {'table', 'detailed', 'compact'}));
    addParameter(p, 'Delimiter', '\t', @ischar);
    parse(p, ConditionNames, NCells, NTraj, CDFs, MSDfitTraj, outputFilename, varargin{:});
    
    % Validate input dimensions
    n = length(ConditionNames);
    if length(NCells) ~= n || length(NTraj) ~= n || length(CDFs) ~= n || length(MSDfitTraj) ~= n
        error('All input arrays must have the same length');
    end
    
    fprintf('Exporting summary statistics for %d conditions...\n', n);
    
    % Calculate CDF lengths
    CDFLengths = zeros(n, 1);
    for i = 1:n
        if ~isempty(CDFs{i}) && isstruct(CDFs{i}) && isfield(CDFs{i}, 'CDFOfJumps')
            CDFLengths(i) = length(CDFs{i}.CDFOfJumps);
        elseif ~isempty(CDFs{i}) && isnumeric(CDFs{i})
            CDFLengths(i) = length(CDFs{i});
        else
            CDFLengths(i) = 0;
            fprintf('Warning: CDF data for condition %d (%s) is empty or invalid\n', i, ConditionNames{i});
        end
    end
    
    % Calculate MSDfitTraj statistics per condition
    MSDfitNvals   = zeros(n, 1);  % number of D estimates
    MSDfitMean    = zeros(n, 1);
    MSDfitMedian  = zeros(n, 1);
    MSDfitStd     = zeros(n, 1);
    MSDfitMin     = zeros(n, 1);
    MSDfitMax     = zeros(n, 1);
    for i = 1:n
        if ~isempty(MSDfitTraj{i}) && isnumeric(MSDfitTraj{i})
            vals = MSDfitTraj{i}(:);
            vals = vals(~isnan(vals));  % exclude NaNs
            MSDfitNvals(i)  = length(vals);
            if MSDfitNvals(i) > 0
                MSDfitMean(i)   = mean(vals);
                MSDfitMedian(i) = median(vals);
                MSDfitStd(i)    = std(vals);
                MSDfitMin(i)    = min(vals);
                MSDfitMax(i)    = max(vals);
            end
        else
            MSDfitNvals(i) = 0;
            fprintf('Warning: MSDfitTraj data for condition %d (%s) is empty or invalid\n', i, ConditionNames{i});
        end
    end
    
    % Pack MSD stats into a struct for cleaner passing
    MSDstats.Nvals  = MSDfitNvals;
    MSDstats.Mean   = MSDfitMean;
    MSDstats.Median = MSDfitMedian;
    MSDstats.Std    = MSDfitStd;
    MSDstats.Min    = MSDfitMin;
    MSDstats.Max    = MSDfitMax;
    
    % Export based on format
    switch p.Results.Format
        case 'table'
            writeTableFormat(ConditionNames, NCells, NTraj, CDFLengths, MSDstats, outputFilename, p.Results.Delimiter);
        case 'detailed'
            writeDetailedFormat(ConditionNames, NCells, NTraj, CDFLengths, MSDstats, outputFilename);
        case 'compact'
            writeCompactFormat(ConditionNames, NCells, NTraj, CDFLengths, MSDstats, outputFilename);
    end
    
    fprintf('Summary statistics exported to: %s\n', outputFilename);
end

function writeTableFormat(ConditionNames, NCells, NTraj, CDFLengths, MSDstats, filename, delimiter)
% Write data in tabular format
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    d = delimiter;
    
    % Write header
    fprintf(fid, 'CDF & MSD-fit Summary Statistics\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, 'Total Conditions: %d\n\n', length(ConditionNames));
    
    % Write column headers
    fprintf(fid, ['Condition%sNCells%sNTraj%sCDF_Length%s' ...
                   'MSDfit_N%sMSDfit_Mean%sMSDfit_Median%sMSDfit_Std%sMSDfit_Min%sMSDfit_Max\n'], ...
                   d, d, d, d, d, d, d, d, d);
    fprintf(fid, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n', ...
        repmat('-',1,20), d, repmat('-',1,10), d, repmat('-',1,10), d, repmat('-',1,12), d, ...
        repmat('-',1,10), d, repmat('-',1,12), d, repmat('-',1,12), d, repmat('-',1,12), d, ...
        repmat('-',1,12), d, repmat('-',1,12));
    
    % Write data rows
    for i = 1:length(ConditionNames)
        fprintf(fid, '%s%s%d%s%d%s%d%s%d%s%.4e%s%.4e%s%.4e%s%.4e%s%.4e\n', ...
            ConditionNames{i}, d, ...
            NCells(i), d, ...
            NTraj(i), d, ...
            CDFLengths(i), d, ...
            MSDstats.Nvals(i), d, ...
            MSDstats.Mean(i), d, ...
            MSDstats.Median(i), d, ...
            MSDstats.Std(i), d, ...
            MSDstats.Min(i), d, ...
            MSDstats.Max(i));
    end
    
    % Write summary statistics
    fprintf(fid, '\n%s\n', repmat('=', 1, 80));
    fprintf(fid, 'SUMMARY STATISTICS\n');
    fprintf(fid, '%s\n', repmat('=', 1, 80));
    fprintf(fid, 'Total Cells:%s%d\n', d, sum(NCells));
    fprintf(fid, 'Total Trajectories:%s%d\n', d, sum(NTraj));
    fprintf(fid, 'Total CDF Points:%s%d\n', d, sum(CDFLengths));
    fprintf(fid, 'Total MSDfit D Estimates:%s%d\n', d, sum(MSDstats.Nvals));
    fprintf(fid, 'Average Cells per Condition:%s%.1f\n', d, mean(NCells));
    fprintf(fid, 'Average Trajectories per Condition:%s%.1f\n', d, mean(NTraj));
    fprintf(fid, 'Average CDF Length:%s%.1f\n', d, mean(CDFLengths));
    fprintf(fid, 'Average MSDfit D Estimates per Condition:%s%.1f\n', d, mean(MSDstats.Nvals));
    fprintf(fid, 'Grand Mean MSDfit D (um^2/s):%s%.4e\n', d, mean(MSDstats.Mean(MSDstats.Nvals > 0)));
    fprintf(fid, 'Grand Median MSDfit D (um^2/s):%s%.4e\n', d, mean(MSDstats.Median(MSDstats.Nvals > 0)));
    
    fclose(fid);
end

function writeDetailedFormat(ConditionNames, NCells, NTraj, CDFLengths, MSDstats, filename)
% Write data in detailed format with descriptions
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    % Write header
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '              CDF & MSD-FIT DIFFUSION SUMMARY STATISTICS\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, 'Total Conditions Analyzed: %d\n', length(ConditionNames));
    fprintf(fid, '\nMSDfitTraj: trajectory-wise diffusion coefficient derived from MSD fit\n\n');
    
    % Write detailed information for each condition
    for i = 1:length(ConditionNames)
        fprintf(fid, 'CONDITION %d: %s\n', i, ConditionNames{i});
        fprintf(fid, '  Number of Cells: %d\n', NCells(i));
        fprintf(fid, '  Number of Trajectories: %d\n', NTraj(i));
        fprintf(fid, '  CDF Data Points: %d\n', CDFLengths(i));
        
        % Calculate derived metrics
        if NCells(i) > 0
            trajPerCell = NTraj(i) / NCells(i);
            fprintf(fid, '  Trajectories per Cell: %.2f\n', trajPerCell);
        end
        
        if NTraj(i) > 0 && CDFLengths(i) > 0
            pointsPerTraj = CDFLengths(i) / NTraj(i);
            fprintf(fid, '  CDF Points per Trajectory: %.1f\n', pointsPerTraj);
        end
        
        % MSDfitTraj statistics
        fprintf(fid, '  --- MSD-fit Diffusion Coefficients ---\n');
        fprintf(fid, '  Number of D Estimates: %d\n', MSDstats.Nvals(i));
        if MSDstats.Nvals(i) > 0
            fprintf(fid, '  Mean D (um^2/s):   %.4e\n', MSDstats.Mean(i));
            fprintf(fid, '  Median D (um^2/s): %.4e\n', MSDstats.Median(i));
            fprintf(fid, '  Std D (um^2/s):    %.4e\n', MSDstats.Std(i));
            fprintf(fid, '  Min D (um^2/s):    %.4e\n', MSDstats.Min(i));
            fprintf(fid, '  Max D (um^2/s):    %.4e\n', MSDstats.Max(i));
        else
            fprintf(fid, '  (no valid D estimates)\n');
        end
        
        fprintf(fid, '\n');
    end
    
    % Write overall summary
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                              OVERALL SUMMARY\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Total Cells across all conditions: %d\n', sum(NCells));
    fprintf(fid, 'Total Trajectories across all conditions: %d\n', sum(NTraj));
    fprintf(fid, 'Total CDF Data Points: %d\n', sum(CDFLengths));
    fprintf(fid, 'Total MSDfit D Estimates: %d\n', sum(MSDstats.Nvals));
    
    fprintf(fid, '\nAverage per Condition:\n');
    fprintf(fid, '  Cells: %.1f (+/-%.1f)\n', mean(NCells), std(NCells));
    fprintf(fid, '  Trajectories: %.1f (+/-%.1f)\n', mean(NTraj), std(NTraj));
    fprintf(fid, '  CDF Points: %.1f (+/-%.1f)\n', mean(CDFLengths), std(CDFLengths));
    fprintf(fid, '  MSDfit D Estimates: %.1f (+/-%.1f)\n', mean(MSDstats.Nvals), std(MSDstats.Nvals));
    
    % MSDfit D summary across conditions (only conditions with data)
    validIdx = MSDstats.Nvals > 0;
    if any(validIdx)
        fprintf(fid, '\nMSD-fit Diffusion Coefficient Summary (across conditions with data):\n');
        fprintf(fid, '  Grand Mean D (um^2/s):   %.4e (+/-%.4e)\n', mean(MSDstats.Mean(validIdx)), std(MSDstats.Mean(validIdx)));
        fprintf(fid, '  Grand Median D (um^2/s): %.4e (+/-%.4e)\n', mean(MSDstats.Median(validIdx)), std(MSDstats.Median(validIdx)));
        fprintf(fid, '  D Mean Range:   %.4e - %.4e\n', min(MSDstats.Mean(validIdx)), max(MSDstats.Mean(validIdx)));
        fprintf(fid, '  D Median Range: %.4e - %.4e\n', min(MSDstats.Median(validIdx)), max(MSDstats.Median(validIdx)));
    end
    
    % Additional statistics
    fprintf(fid, '\nRange Statistics:\n');
    fprintf(fid, '  Cells: %d - %d\n', min(NCells), max(NCells));
    fprintf(fid, '  Trajectories: %d - %d\n', min(NTraj), max(NTraj));
    fprintf(fid, '  CDF Points: %d - %d\n', min(CDFLengths), max(CDFLengths));
    fprintf(fid, '  MSDfit D Estimates: %d - %d\n', min(MSDstats.Nvals), max(MSDstats.Nvals));
    
    fprintf(fid, '\n================================================================================\n');
    
    fclose(fid);
end

function writeCompactFormat(ConditionNames, NCells, NTraj, CDFLengths, MSDstats, filename)
% Write data in compact format
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    % Write compact header
    fprintf(fid, 'CDF & MSDfit Summary (%s) - %d conditions\n', datestr(now, 'yyyy-mm-dd'), length(ConditionNames));
    
    % Write data in compact format
    for i = 1:length(ConditionNames)
        if MSDstats.Nvals(i) > 0
            fprintf(fid, '%s: %d cells, %d traj, %d CDF pts, %d D_est (mean=%.4e, med=%.4e)\n', ...
                ConditionNames{i}, NCells(i), NTraj(i), CDFLengths(i), ...
                MSDstats.Nvals(i), MSDstats.Mean(i), MSDstats.Median(i));
        else
            fprintf(fid, '%s: %d cells, %d traj, %d CDF pts, 0 D_est\n', ...
                ConditionNames{i}, NCells(i), NTraj(i), CDFLengths(i));
        end
    end
    
    % Write totals
    fprintf(fid, 'TOTAL: %d cells, %d traj, %d CDF pts, %d D_est\n', ...
        sum(NCells), sum(NTraj), sum(CDFLengths), sum(MSDstats.Nvals));
    
    fclose(fid);
end

function validateAndDisplaySummary(ConditionNames, NCells, NTraj, CDFs, MSDfitTraj)
% VALIDATEANDDISPLAYSUMMARY Validate inputs and display summary to console
%
% Usage:
%   validateAndDisplaySummary(ConditionNames, NCells, NTraj, CDFs, MSDfitTraj)

    fprintf('\n=== CDF & MSD-fit Data Validation and Summary ===\n');
    
    % Check dimensions
    n = length(ConditionNames);
    if length(NCells) ~= n
        warning('NCells length (%d) does not match ConditionNames length (%d)', length(NCells), n);
    end
    if length(NTraj) ~= n
        warning('NTraj length (%d) does not match ConditionNames length (%d)', length(NTraj), n);
    end
    if length(CDFs) ~= n
        warning('CDFs length (%d) does not match ConditionNames length (%d)', length(CDFs), n);
    end
    if length(MSDfitTraj) ~= n
        warning('MSDfitTraj length (%d) does not match ConditionNames length (%d)', length(MSDfitTraj), n);
    end
    
    % Display summary
    fprintf('Total conditions: %d\n', n);
    fprintf('\nCondition Summary:\n');
    nMin = min([n, length(NCells), length(NTraj), length(CDFs), length(MSDfitTraj)]);
    for i = 1:nMin
        % Calculate CDF length
        if ~isempty(CDFs{i}) && isstruct(CDFs{i}) && isfield(CDFs{i}, 'CDFOfJumps')
            cdfLen = length(CDFs{i}.CDFOfJumps);
        elseif ~isempty(CDFs{i}) && isnumeric(CDFs{i})
            cdfLen = length(CDFs{i});
        else
            cdfLen = 0;
        end
        
        % Calculate MSDfitTraj count
        if ~isempty(MSDfitTraj{i}) && isnumeric(MSDfitTraj{i})
            msdN = sum(~isnan(MSDfitTraj{i}(:)));
        else
            msdN = 0;
        end
        
        fprintf('  %d. %s: %d cells, %d traj, %d CDF pts, %d D_est\n', ...
            i, ConditionNames{i}, NCells(i), NTraj(i), cdfLen, msdN);
    end
    
    if nMin == n
        fprintf('\nTotals: %d cells, %d trajectories\n', sum(NCells(1:n)), sum(NTraj(1:n)));
    end
    
    fprintf('=================================================\n\n');
end