% This function accepts an SMF structure as an input and writes the SMF
% parameters to a .csv file with one row per parameter. Nested struct
% fields are flattened into dotted names (e.g. Thresholding.AutoThreshLogL).
% Columns: Parameter, Value, Class, Size
%
% Short numeric vectors/matrices are written as a single quoted cell such as
% "[1 2 3]" (matrix rows separated by ";"). Arrays with more than 20
% elements are summarized (size and min/max) instead of listed.
%
% Created by Ian Adams, 2025.

% Example usage
%
%  load a saved SMF.mat file
% >> IMA_exportStructToCSV(SMF, "C:\path\to\your\data\SMF_params.csv")
%
% scripted call, need to execute in a different .m file
%
% SMFdir = "C:\Users\imadams\Documents\GroundZero\SMF_optimization";
% filename = "SMF_DNA-PAINT_Atto665_04";
% fileIn = fullfile(SMFdir, strcat(filename, ".mat"));
% fileOut = fullfile(SMFdir, strcat(filename, ".csv"));
%
% SMF = load(fileIn);
% IMA_exportStructToCSV(SMF, fileOut)


function IMA_exportStructToCSV(SMF, filename)
    fileID = fopen(filename, 'w');
    if fileID < 0
        error('IMA_exportStructToCSV:fopen', 'Could not open %s for writing.', filename);
    end
    cleanup = onCleanup(@() fclose(fileID));

    fprintf(fileID, 'Parameter,Value,Class,Size\n');
    writeStruct(fileID, SMF, '');
end

function writeStruct(fileID, S, prefix)
    % Recursively write all fields of struct S with dotted parameter names.
    maxElements = 20;  % arrays larger than this are summarized, not listed
    fields = fieldnames(S);
    for i = 1:numel(fields)
        name = fields{i};
        value = S.(name);
        fullName = [prefix name];

        if isstruct(value) && isscalar(value)
            writeStruct(fileID, value, [fullName '.']);
            continue
        end

        if islogical(value) && isscalar(value)
            if value, str = 'true'; else, str = 'false'; end
        elseif islogical(value) || isnumeric(value)
            if isempty(value)
                str = '[]';
            elseif isscalar(value)
                str = sprintf('%.6g', value);
            elseif numel(value) > maxElements
                str = sprintf('[array, %d elements; min %.6g, max %.6g] (not exported)', ...
                    numel(value), double(min(value(:))), double(max(value(:))));
            elseif isvector(value)
                str = ['[' strtrim(sprintf('%.6g ', value)) ']'];
            elseif ismatrix(value)
                rows = cell(1, size(value, 1));
                for r = 1:size(value, 1)
                    rows{r} = strtrim(sprintf('%.6g ', value(r, :)));
                end
                str = ['[' strjoin(rows, '; ') ']'];
            else
                continue
            end
        elseif ischar(value)
            str = value;
        elseif isstring(value) && isscalar(value)
            str = char(value);
        else
            continue  % cells, struct arrays, etc. are skipped
        end

        sz = sprintf('%dx', size(value));
        fprintf(fileID, '%s,%s,%s,%s\n', csvEscape(fullName), csvEscape(str), ...
            class(value), sz(1:end-1));
    end
end

function out = csvEscape(str)
    % Quote a field if it contains a comma, quote or newline.
    str = strrep(str, newline, ' ');
    if any(str == ',') || any(str == '"')
        out = ['"' strrep(str, '"', '""') '"'];
    else
        out = str;
    end
end
