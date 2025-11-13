function IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, outputFilename, varargin)
% EXPORTCDFSUMMARYSTATS Export CDF summary statistics to text file
%
% Usage:
%   IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, 'summary.txt')
%   IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, 'summary.txt', 'Format', 'table')
%   IMA_exportD_est_Nvals(ConditionNames, NCells, NTraj, CDFs, 'summary.txt', 'Format', 'detailed')
%
% Inputs:
%   ConditionNames  - 1×n cell array of condition names
%   NCells          - n×1 vector of number of cells per condition
%   NTraj           - n×1 vector of number of trajectories per condition
%   CDFs            - n×1 cell array of CDF data structures
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
    addRequired(p, 'outputFilename', @ischar);
    addParameter(p, 'Format', 'table', @(x) ismember(x, {'table', 'detailed', 'compact'}));
    addParameter(p, 'Delimiter', '\t', @ischar);
    parse(p, ConditionNames, NCells, NTraj, CDFs, outputFilename, varargin{:});
    
    % Validate input dimensions
    n = length(ConditionNames);
    if length(NCells) ~= n || length(NTraj) ~= n || length(CDFs) ~= n
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
    
    % Export based on format
    switch p.Results.Format
        case 'table'
            writeTableFormat(ConditionNames, NCells, NTraj, CDFLengths, outputFilename, p.Results.Delimiter);
        case 'detailed'
            writeDetailedFormat(ConditionNames, NCells, NTraj, CDFLengths, outputFilename);
        case 'compact'
            writeCompactFormat(ConditionNames, NCells, NTraj, CDFLengths, outputFilename);
    end
    
    fprintf('Summary statistics exported to: %s\n', outputFilename);
end

function writeTableFormat(ConditionNames, NCells, NTraj, CDFLengths, filename, delimiter)
% Write data in tabular format
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    % Write header
    fprintf(fid, 'CDF Summary Statistics\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, 'Total Conditions: %d\n\n', length(ConditionNames));
    
    % Write column headers
    fprintf(fid, 'Condition%sNCells%sNTraj%sCDF_Length\n', delimiter, delimiter, delimiter);
    fprintf(fid, '%s%s%s%s%s\n', repmat('-', 1, 20), delimiter, repmat('-', 1, 10), delimiter, repmat('-', 1, 10), delimiter, repmat('-', 1, 15));
    
    % Write data rows
    for i = 1:length(ConditionNames)
        fprintf(fid, '%s%s%d%s%d%s%d\n', ...
            ConditionNames{i}, delimiter, ...
            NCells(i), delimiter, ...
            NTraj(i), delimiter, ...
            CDFLengths(i));
    end
    
    % Write summary statistics
    fprintf(fid, '\n%s\n', repmat('=', 1, 60));
    fprintf(fid, 'SUMMARY STATISTICS\n');
    fprintf(fid, '%s\n', repmat('=', 1, 60));
    fprintf(fid, 'Total Cells:%s%d\n', delimiter, sum(NCells));
    fprintf(fid, 'Total Trajectories:%s%d\n', delimiter, sum(NTraj));
    fprintf(fid, 'Total CDF Points:%s%d\n', delimiter, sum(CDFLengths));
    fprintf(fid, 'Average Cells per Condition:%s%.1f\n', delimiter, mean(NCells));
    fprintf(fid, 'Average Trajectories per Condition:%s%.1f\n', delimiter, mean(NTraj));
    fprintf(fid, 'Average CDF Length:%s%.1f\n', delimiter, mean(CDFLengths));
    
    fclose(fid);
end

function writeDetailedFormat(ConditionNames, NCells, NTraj, CDFLengths, filename)
% Write data in detailed format with descriptions
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    % Write header
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                           CDF SUMMARY STATISTICS\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, 'Total Conditions Analyzed: %d\n\n', length(ConditionNames));
    
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
        
        fprintf(fid, '\n');
    end
    
    % Write overall summary
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                              OVERALL SUMMARY\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Total Cells across all conditions: %d\n', sum(NCells));
    fprintf(fid, 'Total Trajectories across all conditions: %d\n', sum(NTraj));
    fprintf(fid, 'Total CDF Data Points: %d\n', sum(CDFLengths));
    fprintf(fid, '\nAverage per Condition:\n');
    fprintf(fid, '  Cells: %.1f (±%.1f)\n', mean(NCells), std(NCells));
    fprintf(fid, '  Trajectories: %.1f (±%.1f)\n', mean(NTraj), std(NTraj));
    fprintf(fid, '  CDF Points: %.1f (±%.1f)\n', mean(CDFLengths), std(CDFLengths));
    
    % Additional statistics
    fprintf(fid, '\nRange Statistics:\n');
    fprintf(fid, '  Cells: %d - %d\n', min(NCells), max(NCells));
    fprintf(fid, '  Trajectories: %d - %d\n', min(NTraj), max(NTraj));
    fprintf(fid, '  CDF Points: %d - %d\n', min(CDFLengths), max(CDFLengths));
    
    fprintf(fid, '\n================================================================================\n');
    
    fclose(fid);
end

function writeCompactFormat(ConditionNames, NCells, NTraj, CDFLengths, filename)
% Write data in compact format
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    % Write compact header
    fprintf(fid, 'CDF Summary (%s) - %d conditions\n', datestr(now, 'yyyy-mm-dd'), length(ConditionNames));
    
    % Write data in compact format
    for i = 1:length(ConditionNames)
        fprintf(fid, '%s: %d cells, %d traj, %d CDF pts\n', ...
            ConditionNames{i}, NCells(i), NTraj(i), CDFLengths(i));
    end
    
    % Write totals
    fprintf(fid, 'TOTAL: %d cells, %d traj, %d CDF pts\n', ...
        sum(NCells), sum(NTraj), sum(CDFLengths));
    
    fclose(fid);
end

function validateAndDisplaySummary(ConditionNames, NCells, NTraj, CDFs)
% VALIDATEANDDISPLAYSUMMARY Validate inputs and display summary to console
%
% Usage:
%   validateAndDisplaySummary(ConditionNames, NCells, NTraj, CDFs)

    fprintf('\n=== CDF Data Validation and Summary ===\n');
    
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
    
    % Display summary
    fprintf('Total conditions: %d\n', n);
    fprintf('\nCondition Summary:\n');
    for i = 1:min(n, length(NCells), length(NTraj), length(CDFs))
        % Calculate CDF length
        if ~isempty(CDFs{i}) && isstruct(CDFs{i}) && isfield(CDFs{i}, 'CDFOfJumps')
            cdfLen = length(CDFs{i}.CDFOfJumps);
        elseif ~isempty(CDFs{i}) && isnumeric(CDFs{i})
            cdfLen = length(CDFs{i});
        else
            cdfLen = 0;
        end
        
        fprintf('  %d. %s: %d cells, %d traj, %d CDF points\n', ...
            i, ConditionNames{i}, NCells(i), NTraj(i), cdfLen);
    end
    
    if n <= length(NCells) && n <= length(NTraj) && n <= length(CDFs)
        fprintf('\nTotals: %d cells, %d trajectories\n', sum(NCells(1:n)), sum(NTraj(1:n)));
    end
    
    fprintf('=====================================\n\n');
end