% This function accepts an SMF structure as an input and writes the SMF
% parameters to a text document. This is helpful to be able to save some
% notes about your analyses outside of the SMF.m file that is typically
% saved with the analysis parameters for SMLM work. 

% Created by Ian Adams, 2025. 

% Example usage
%
%  load  a saved SMF.mat file
%  In the Matlab command window, call the function. The second argument is
%  the full path for the desired save location with ".txt" extension.
% >> IMA_exportStructToText(SMF, "C:\path\to\your\data\SMF_text.txt")
% 
%
% scripted call, need to execute in a different .m file
%
% SMFdir = "C:\Users\imadams\Documents\GroundZero\SMF_optimization";
% filename = "SMF_DNA-PAINT_Atto665_04";
% fileIn = fullfile(SMFdir, strcat(filename, ".mat"));
% fileOut = fullfile(SMFdir, strcat(filename, ".txt"));
% 
% SMF = load(fileIn);
% IMA_exportStructToText(SMF, fileOut)


function exportStructToText(SMF, filename)
    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Write header
    fprintf(fileID, 'SMF Parameters Export\n');
    fprintf(fileID, '===================\n\n');

    writeStruct(fileID, SMF, 0);

    % Close the file
    fclose(fileID);
end

function writeStruct(fileID, S, depth)
    % Recursively write all fields of struct S, indenting by depth.
    maxElements = 20;  % arrays larger than this are summarized, not listed
    indent = repmat('\t', 1, depth);
    fields = fieldnames(S);
    for i = 1:numel(fields)
        name = fields{i};
        value = S.(name);

        if isstruct(value) && isscalar(value)
            fprintf(fileID, [indent '%s (nested struct):\n'], name);
            writeStruct(fileID, value, depth + 1);
        elseif islogical(value) && isscalar(value)
            % Must be checked alongside numeric: isnumeric(true) is false
            if value
                fprintf(fileID, [indent '%s: true\n'], name);
            else
                fprintf(fileID, [indent '%s: false\n'], name);
            end
        elseif islogical(value) || isnumeric(value)
            if isempty(value)
                fprintf(fileID, [indent '%s: []\n'], name);
            elseif isscalar(value)
                fprintf(fileID, [indent '%s: %.6f\n'], name, value);
            elseif numel(value) > maxElements
                % Large data array: write a summary instead of every value
                sz = sprintf('%dx', size(value));
                fprintf(fileID, [indent '%s: [%s %s array, %d elements; min %.6g, max %.6g] (not exported)\n'], ...
                    name, sz(1:end-1), class(value), numel(value), ...
                    double(min(value(:))), double(max(value(:))));
            elseif isvector(value)
                fprintf(fileID, [indent '%s: ['], name);
                fprintf(fileID, '%.6f ', value);
                fprintf(fileID, ']\n');
            elseif ismatrix(value)
                fprintf(fileID, [indent '%s:\n'], name);
                for row = 1:size(value, 1)
                    fprintf(fileID, [indent '\t']);
                    fprintf(fileID, '%.6f ', value(row, :));
                    fprintf(fileID, '\n');
                end
            end
        elseif ischar(value)
            fprintf(fileID, [indent '%s: %s\n'], name, value);
        elseif isstring(value) && isscalar(value)
            fprintf(fileID, [indent '%s: %s\n'], name, char(value));
        end
    end
end
