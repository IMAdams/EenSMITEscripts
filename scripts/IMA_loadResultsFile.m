function [TR1, SMF1, TR2, SMF2] = IMA_loadResultsFile(filePath1, filePath2)
    % Load the first .mat file
    data1 = load(filePath1);
    TR1 = data1.TR;
    SMF1 = data1.SMF;
    
    % Load the second .mat file
    data2 = load(filePath2);
    TR2 = data2.TR;
    SMF2 = data2.SMF;
    
    % Error checking
    if ~exist('TR1', 'var') || ~exist('SMF1', 'var') || ...
       ~exist('TR2', 'var') || ~exist('SMF2', 'var')
        error('One or both results files do not contain TR and SMF structures');
    end
end