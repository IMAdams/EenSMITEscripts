function IMA_extractCDFJumpsToh5(CDFEnsemble, filename, varargin)
% EXPORTCDFTOHDF5 Export CDFEnsemble data to HDF5 format for memory efficiency
%
% Usage:
%   exportCDFToHDF5(CDFEnsemble, 'data.h5')
%   exportCDFToHDF5(CDFEnsemble, 'data.h5', 'DatasetNames', {'Exp1', 'Exp2'})
%   exportCDFToHDF5(CDFEnsemble, 'data.h5', 'FieldsToExport', {'CDFOfJumps', 'SortedSquaredDisp'})
%
% Inputs:
%   CDFEnsemble     - Cell array of structs with CDF data
%   filename        - String, name of output HDF5 file
%
% Optional Parameters:
%   'DatasetNames'   - Cell array of dataset identifiers
%   'FieldsToExport' - Cell array of field names to export (default: all CDF fields)
%   'Compression'    - HDF5 compression level 0-9 (default: 6)

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'CDFEnsemble');
    addRequired(p, 'filename', @ischar);
    addParameter(p, 'DatasetNames', {}, @iscell);
    addParameter(p, 'FieldsToExport', {'CDFOfJumps', 'SortedSquaredDisp', 'SquaredDisplacement', 'FrameLagsAll', 'LocVarianceSum', 'CellID'}, @iscell);
    addParameter(p, 'Compression', 6, @(x) isnumeric(x) && x >= 0 && x <= 9);
    parse(p, CDFEnsemble, filename, varargin{:});
    
    % Ensure CDFEnsemble is a cell array
    if ~iscell(CDFEnsemble)
        CDFEnsemble = {CDFEnsemble};
    end
    
    % Delete existing file if it exists
    if exist(filename, 'file')
        delete(filename);
        fprintf('Deleted existing file: %s\n', filename);
    end
    
    fprintf('Exporting %d datasets to HDF5 format...\n', length(CDFEnsemble));
    
    % Export each dataset
    for i = 1:length(CDFEnsemble)
        data = CDFEnsemble{i};
        
        % Dataset identifier
        if ~isempty(p.Results.DatasetNames) && length(p.Results.DatasetNames) >= i
            datasetName = p.Results.DatasetNames{i};
        else
            datasetName = sprintf('Dataset_%d', i);
        end
        
        fprintf('  Processing %s...\n', datasetName);
        
        % Create group for this dataset
        groupName = ['/' datasetName];
        
        % Export each requested field
        for j = 1:length(p.Results.FieldsToExport)
            fieldName = p.Results.FieldsToExport{j};
            
            if isfield(data, fieldName) && ~isempty(data.(fieldName))
                fieldData = data.(fieldName);
                datasetPath = [groupName '/' fieldName];
                
                % Handle different data types
                if ischar(fieldData) || isstring(fieldData)
                    % Handle string/char data
                    try
                        fieldData = char(fieldData); % Convert string to char if needed
                        h5create(filename, datasetPath, size(fieldData), ...
                            'Datatype', 'string');
                        h5write(filename, datasetPath, fieldData);
                        
                        % Add string-specific attributes
                        h5writeatt(filename, datasetPath, 'description', ...
                            sprintf('%s string data from %s', fieldName, datasetName));
                        h5writeatt(filename, datasetPath, 'data_type', 'string');
                        h5writeatt(filename, datasetPath, 'string_length', length(fieldData));
                        
                    catch ME
                        fprintf('    Warning: Could not export string field %s - %s\n', fieldName, ME.message);
                    end
                    
                elseif iscell(fieldData)
                    % Handle cell arrays (e.g., cell array of strings)
                    try
                        if all(cellfun(@(x) ischar(x) || isstring(x), fieldData))
                            % Cell array of strings - convert to string array
                            stringArray = string(fieldData);
                            h5create(filename, datasetPath, size(stringArray), ...
                                'Datatype', 'string');
                            h5write(filename, datasetPath, stringArray);
                            
                            % Add cell array attributes
                            h5writeatt(filename, datasetPath, 'description', ...
                                sprintf('%s cell array of strings from %s', fieldName, datasetName));
                            h5writeatt(filename, datasetPath, 'data_type', 'cell_array_of_strings');
                            h5writeatt(filename, datasetPath, 'array_length', length(fieldData));
                        else
                            fprintf('    Warning: %s contains non-string cell array data - skipping\n', fieldName);
                        end
                    catch ME
                        fprintf('    Warning: Could not export cell array field %s - %s\n', fieldName, ME.message);
                    end
                    
                elseif isnumeric(fieldData)
                    % Handle numeric data (original functionality)
                    try
                        h5create(filename, datasetPath, size(fieldData), ...
                            'Datatype', 'double', ...
                            'ChunkSize', min(size(fieldData), [10000, 1]), ...
                            'Deflate', p.Results.Compression);
                        h5write(filename, datasetPath, fieldData);
                        
                        % Add numeric-specific attributes
                        h5writeatt(filename, datasetPath, 'description', ...
                            sprintf('%s numeric data from %s', fieldName, datasetName));
                        h5writeatt(filename, datasetPath, 'length', length(fieldData));
                        h5writeatt(filename, datasetPath, 'data_type', 'numeric');
                        
                    catch ME
                        fprintf('    Warning: Could not export numeric field %s - %s\n', fieldName, ME.message);
                    end
                else
                    fprintf('    Warning: %s has unsupported data type (%s) - skipping\n', fieldName, class(fieldData));
                end
            else
                fprintf('    Warning: Field %s not found or empty in %s\n', fieldName, datasetName);
            end
        end
        
        % Add dataset-level metadata
        try
            h5writeatt(filename, groupName, 'dataset_index', i);
            h5writeatt(filename, groupName, 'export_timestamp', datestr(now));
            h5writeatt(filename, groupName, 'matlab_version', version);
        catch
            % Attributes are optional, continue if they fail
        end
    end
    
    % Add file-level metadata
    try
        h5writeatt(filename, '/', 'total_datasets', length(CDFEnsemble));
        h5writeatt(filename, '/', 'export_date', datestr(now));
        h5writeatt(filename, '/', 'exported_fields', strjoin(p.Results.FieldsToExport, ', '));
        h5writeatt(filename, '/', 'compression_level', p.Results.Compression);
    catch
        % Attributes are optional
    end
    
    fprintf('Export complete: %s\n', filename);
    displayHDF5Info(filename);
end

function displayHDF5Info(filename)
% Display information about the created HDF5 file
    try
        info = h5info(filename);
        fileInfo = dir(filename);
        
        fprintf('\nHDF5 File Information:\n');
        fprintf('  File: %s\n', filename);
        fprintf('  Size: %.2f MB\n', fileInfo.bytes / 1024^2);
        fprintf('  Groups: %d\n', length(info.Groups));
        
        for i = 1:length(info.Groups)
            group = info.Groups(i);
            fprintf('    %s: %d datasets\n', group.Name, length(group.Datasets));
        end
        
    catch ME
        fprintf('Could not display file info: %s\n', ME.message);
    end
end

function data = readCDFFromHDF5(filename, varargin)
% READCDFFROMHDF5 Read CDF data back from HDF5 file
%
% Usage:
%   data = readCDFFromHDF5('data.h5')  % Read all datasets and fields
%   data = readCDFFromHDF5('data.h5', 'Dataset', 'Dataset_1')  % Read specific dataset
%   data = readCDFFromHDF5('data.h5', 'Fields', {'CDFOfJumps'})  % Read specific fields
%
% Returns:
%   data - Struct or cell array of structs with the CDF data

    p = inputParser;
    addRequired(p, 'filename', @ischar);
    addParameter(p, 'Dataset', '', @ischar);  % Specific dataset name
    addParameter(p, 'Fields', {}, @iscell);   % Specific fields to read
    parse(p, filename, varargin{:});
    
    if ~exist(filename, 'file')
        error('HDF5 file does not exist: %s', filename);
    end
    
    try
        info = h5info(filename);
    catch ME
        error('Cannot read HDF5 file info: %s', ME.message);
    end
    
    % Determine which datasets to read
    if ~isempty(p.Results.Dataset)
        % Read specific dataset
        targetGroups = {p.Results.Dataset};
        groupsToRead = {};
        for i = 1:length(info.Groups)
            groupName = info.Groups(i).Name(2:end); % Remove leading '/'
            if strcmp(groupName, p.Results.Dataset)
                groupsToRead{end+1} = info.Groups(i);
                break;
            end
        end
        if isempty(groupsToRead)
            error('Dataset not found: %s', p.Results.Dataset);
        end
    else
        % Read all datasets
        groupsToRead = info.Groups;
    end
    
    fprintf('Reading %d datasets from HDF5 file...\n', length(groupsToRead));
    
    % Read data
    if length(groupsToRead) == 1 && ~isempty(p.Results.Dataset)
        % Return single struct for specific dataset request
        data = readGroupData(filename, groupsToRead{1}, p.Results.Fields);
    else
        % Return cell array of structs
        data = {};
        for i = 1:length(groupsToRead)
            group = groupsToRead{i};
            groupName = group.Name(2:end); % Remove leading '/'
            fprintf('  Reading %s...\n', groupName);
            data{end+1} = readGroupData(filename, group, p.Results.Fields);
        end
    end
    
    fprintf('Read complete.\n');
end

function groupData = readGroupData(filename, group, fieldsToRead)
% Read data from a specific HDF5 group
    groupData = struct();
    
    % Determine which datasets to read
    if isempty(fieldsToRead)
        datasetsToRead = group.Datasets;
    else
        datasetsToRead = {};
        for i = 1:length(group.Datasets)
            if ismember(group.Datasets(i).Name, fieldsToRead)
                datasetsToRead{end+1} = group.Datasets(i);
            end
        end
    end
    
    % Read each dataset
    for i = 1:length(datasetsToRead)
        dataset = datasetsToRead{i};
        datasetPath = [group.Name '/' dataset.Name];
        
        try
            fieldData = h5read(filename, datasetPath);
            groupData.(dataset.Name) = fieldData;
        catch ME
            fprintf('    Warning: Could not read %s - %s\n', dataset.Name, ME.message);
        end
    end
end

function exportHDF5ToCSV(hdf5Filename, csvFilename, varargin)
% EXPORTHDF5TOCSV Convert HDF5 CDF data back to CSV format
%
% Usage:
%   exportHDF5ToCSV('data.h5', 'output.csv')
%   exportHDF5ToCSV('data.h5', 'output.csv', 'Dataset', 'Dataset_1')

    p = inputParser;
    addRequired(p, 'hdf5Filename', @ischar);
    addRequired(p, 'csvFilename', @ischar);
    addParameter(p, 'Dataset', '', @ischar);
    parse(p, hdf5Filename, csvFilename, varargin{:});
    
    % Read data from HDF5
    fprintf('Reading HDF5 data...\n');
    if ~isempty(p.Results.Dataset)
        data = readCDFFromHDF5(hdf5Filename, 'Dataset', p.Results.Dataset);
        datasets = {data};
        datasetNames = {p.Results.Dataset};
    else
        datasets = readCDFFromHDF5(hdf5Filename);
        datasetNames = {};
        for i = 1:length(datasets)
            datasetNames{i} = sprintf('Dataset_%d', i);
        end
    end
    
    % Convert to CSV using chunked writing
    fprintf('Converting to CSV...\n');
    fid = fopen(csvFilename, 'w');
    if fid == -1
        error('Cannot open CSV file for writing: %s', csvFilename);
    end
    
    % Write headers
    headers = {'Dataset_ID', 'Point_Index', 'CDF_Value', 'SortedSquaredDisp', 'SquaredDisplacement', 'FrameLag', 'CellID'};
    fprintf(fid, '%s\n', strjoin(headers, ','));
    
    % Write data
    for i = 1:length(datasets)
        data = datasets{i};
        datasetID = datasetNames{i};
        
        if isfield(data, 'CDFOfJumps')
            dataLength = length(data.CDFOfJumps);
            
            % Get CellID for this dataset (single string or cell array)
            cellID = '';
            if isfield(data, 'CellID')
                if ischar(data.CellID) || isstring(data.CellID)
                    cellID = char(data.CellID);
                elseif iscell(data.CellID) && ~isempty(data.CellID)
                    cellID = char(data.CellID{1}); % Use first cell ID if it's a cell array
                end
            end
            
            % Write in chunks to avoid memory issues
            chunkSize = 10000;
            for startIdx = 1:chunkSize:dataLength
                endIdx = min(startIdx + chunkSize - 1, dataLength);
                
                for j = startIdx:endIdx
                    fprintf(fid, '%s,%d,%.6g,%.6g,%.6g,%d,%s\n', ...
                        datasetID, j, ...
                        getFieldValue(data, 'CDFOfJumps', j), ...
                        getFieldValue(data, 'SortedSquaredDisp', j), ...
                        getFieldValue(data, 'SquaredDisplacement', j), ...
                        round(getFieldValue(data, 'FrameLagsAll', j)), ...
                        cellID);
                end
            end
        end
    end
    
    fclose(fid);
    fprintf('CSV export complete: %s\n', csvFilename);
end

function value = getFieldValue(data, fieldName, index)
% Safely get field value with fallback to NaN
    if isfield(data, fieldName) && length(data.(fieldName)) >= index
        value = data.(fieldName)(index);
    else
        value = NaN;
    end
end