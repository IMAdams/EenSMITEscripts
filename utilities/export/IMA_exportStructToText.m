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
    
    % Get all field names from the struct
    % if isstruct(SMF)
        fields = fieldnames(SMF);
    % else
    %     fields = 1;
    % end
    
    % Write header
    fprintf(fileID, 'SMF Parameters Export\n');
    fprintf(fileID, '===================\n\n');
    
    % Iterate through all fields
    for i = 1:length(fields)
        currentField = fields{i};
        value = SMF.(currentField);
        
        % Handle different data types
        if isnumeric(value)
            if isscalar(value)
                % For single numbers
                fprintf(fileID, '%s: %.6f\n', currentField, value);
            elseif isvector(value)
                % For vectors
                fprintf(fileID, '%s: [', currentField);
                fprintf(fileID, '%.6f ', value);
                fprintf(fileID, ']\n');
            elseif ismatrix(value)
                % For matrices
                fprintf(fileID, '%s:\n', currentField);
                for row = 1:size(value, 1)
                    fprintf(fileID, '\t');
                    fprintf(fileID, '%.6f ', value(row, :));
                    fprintf(fileID, '\n');
                end
            end
        elseif ischar(value)
            % For strings
            fprintf(fileID, '%s: %s\n', currentField, value);
        elseif islogical(value)
            % For boolean values
            if value
                fprintf(fileID, '%s: true\n', currentField);
            else
                fprintf(fileID, '%s: false\n', currentField);
            end
        elseif isstruct(value)
            % For nested structs
            fprintf(fileID, '%s (nested struct):\n', currentField);
            nestedFields = fieldnames(value);
            for j = 1:length(nestedFields)
                nestedValue = value.(nestedFields{j});
                if isnumeric(nestedValue) && isscalar(nestedValue)
                    fprintf(fileID, '\t%s: %.6f\n', nestedFields{j}, nestedValue);
                elseif ischar(nestedValue)
                    fprintf(fileID, '\t%s: %s\n', nestedFields{j}, nestedValue);
                end
            end
        end
    end
    
    % Close the file
    fclose(fileID);
end
