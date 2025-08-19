function IMA_extractDataToCSV(CDFEnsemble, outputFilename, varargin)
% EXTRACTCDFJUMPSTOCSV Extract CDFOfJumps data from CDFEnsemble and save to CSV
%
% Usage:
%   extractCDFJumpsToCSV(CDFEnsemble, 'output.csv')
%   extractCDFJumpsToCSV(CDFEnsemble, 'output.csv', 'DatasetNames', {'Exp1', 'Exp2'})
%
% Inputs:
%   CDFEnsemble     - Cell array of structs with CDFOfJumps field
%   outputFilename  - String, name of output CSV file
%
% Optional Parameters:
%   'DatasetNames'   - Cell array of dataset identifiers

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'CDFEnsemble');
    addRequired(p, 'outputFilename', @ischar);
    addParameter(p, 'DatasetNames', {}, @iscell);
    parse(p, CDFEnsemble, outputFilename, varargin{:});
    
    % Ensure CDFEnsemble is a cell array
    if ~iscell(CDFEnsemble)
        CDFEnsemble = {CDFEnsemble};
    end
    
    fprintf('Extracting CDFOfJumps data from %d datasets...\n', length(CDFEnsemble));
    
    % Extract CDF jumps data
    [extractedData, headers] = extractCDFJumpsData(CDFEnsemble, p.Results);
    
    % Write to CSV
    writeDataToCSV(extractedData, headers, outputFilename);
    fprintf('Data exported to: %s\n', outputFilename);
    fprintf('Exported %d rows and %d columns\n', size(extractedData, 1), size(extractedData, 2));
end

function [extractedData, headers] = extractCDFJumpsData(CDFEnsemble, options)
% Extract raw CDFOfJumps data with related displacement information
    extractedData = {};
    
    % Process each dataset
    for i = 1:length(CDFEnsemble)
        data = CDFEnsemble{i};
        
        % Dataset identifier
        if ~isempty(options.DatasetNames) && length(options.DatasetNames) >= i
            datasetID = options.DatasetNames{i};
        else
            datasetID = sprintf('Dataset_%d', i);
        end
        
        if isfield(data, 'CDFOfJumps')
            cdfJumps = data.CDFOfJumps;
            
            % Also extract related fields if available
            if isfield(data, 'SortedSquaredDisp')
                sortedDisp = data.SortedSquaredDisp;
            else
                sortedDisp = [];
            end
            
            if isfield(data, 'SquaredDisplacement')
                squaredDisp = data.SquaredDisplacement;
            else
                squaredDisp = [];
            end
            
            if isfield(data, 'FrameLagsAll')
                frameLagsAll = data.FrameLagsAll;
            else
                frameLagsAll = [];
            end
            
            % Create one row per CDF point
            for j = 1:length(cdfJumps)
                rowData = {datasetID, j, cdfJumps(j)};
                
                % Add related displacement data if available
                if ~isempty(sortedDisp) && j <= length(sortedDisp)
                    rowData{end+1} = sortedDisp(j);
                else
                    rowData{end+1} = NaN;
                end
                
                if ~isempty(squaredDisp) && j <= length(squaredDisp)
                    rowData{end+1} = squaredDisp(j);
                else
                    rowData{end+1} = NaN;
                end
                
                if ~isempty(frameLagsAll) && j <= length(frameLagsAll)
                    rowData{end+1} = frameLagsAll(j);
                else
                    rowData{end+1} = NaN;
                end
                
                extractedData = [extractedData; rowData];
            end
        else
            fprintf('Warning: Dataset %d does not contain CDFOfJumps field\n', i);
        end
    end
    
    headers = {'Dataset_ID', 'Point_Index', 'CDF_Value', 'SortedSquaredDisp', 'SquaredDisplacement', 'FrameLag'};
end

function writeDataToCSV(data, headers, filename)
% Write cell array data to CSV file
    try
        % Convert data to table for easier CSV writing
        dataTable = cell2table(data, 'VariableNames', headers);
        writetable(dataTable, filename);
    catch ME
        fprintf('Table method failed: %s\nUsing fallback method...\n', ME.message);
        % Fallback method using low-level file operations
        fid = fopen(filename, 'w');
        if fid == -1
            error('Could not open file for writing: %s', filename);
        end
        
        % Write headers
        fprintf(fid, '%s', headers{1});
        for i = 2:length(headers)
            fprintf(fid, ',%s', headers{i});
        end
        fprintf(fid, '\n');
        
        % Write data
        for i = 1:size(data, 1)
            for j = 1:size(data, 2)
                if j > 1
                    fprintf(fid, ',');
                end
                
                if isnumeric(data{i,j})
                    if isnan(data{i,j})
                        fprintf(fid, '');
                    else
                        fprintf(fid, '%.6g', data{i,j});
                    end
                else
                    fprintf(fid, '%s', char(data{i,j}));
                end
            end
            fprintf(fid, '\n');
        end
        
        fclose(fid);
    end
end