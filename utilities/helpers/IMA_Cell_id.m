function out = IMA_Cell_id(filename)
    % Extract the first two characters from a filename
    % Input: filename - string containing the full file path or just filename
    % Output: firstTwoChars - string containing the first two characters
    
    % Extract just the filename if a full path is provided
    [~, name, ext] = fileparts(filename);
    fileName = string(name + ext);
    
    % Get cell id
    if strlength(name) >= 18
        out = extractBefore(name, strlength(name)-16);
    else
        out = name;  % Return the whole string if less than 2 chars
    end
end