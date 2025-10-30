%% Define directories, determine files to be analyzed and separated
%%%IMPORTANT: This script only works on 1 or 2 channel data from the Olympus IX83%%%
%%%When you choose files, only choose appropriate files in groups for the sections below

% Define the directory containing the raw data.
%RawDataDir = 'O:\Cell Path\Lidke Lab\Emmy\Data\SiMPull_TCells\20250613_PD-1peptides-anti-alfa488_EGFR-GFP_PD-1-GFP';
RawDataDir = ...
    'O:\Cell Path\Lidke Lab\IMAdams\Data\HeLa-ALFA-EGFR-KI\HelaALFAEGFR_10242025';
%allfiles = ...
%    selectMultipleFiles(RawDataDir, '.vsi', 'Select vsi Files to Combine');
allfiles = uipickfiles('FilterSpec', RawDataDir, 'REFilter', '\.vsi', ...
                       'Prompt', 'Select vsi Files to Combine')';
% Define a save directory in which the tracking results will be stored.
SaveDir = fullfile(RawDataDir, 'data_dot_mat');
if not(isfolder(SaveDir))
    mkdir(SaveDir);
end

%SUGGESTION Can we add something to decifer if the second channel is 561 or 647?
%Essenitally, if we were able to do this we could ensure the data setjoined
%knew to put the second cahannel int he top left quadrant or bottom left.
%This only matter is we move to actual 3 channel data or to avoid
%confustion with 647 data potentially being labeled with 561 in the future.


% Data formats:Run section based on your data type
%    Section A) 2-channel sequence data, ie Fiducial, gain, background
%    taken 20 frames
%    Section B) 2-channel data, single frame- ie TIRF 488/647 channels for
%    SiMPull (each vsi contains both channels of data)
%    Section C) 2-channel sequence data taken sequentially, ie 488/561
%    SiMPull data where each channel has its own vsi file
%    Section D) 1-channel snapshot data
%    

%% SECTION A) 2 Channel Sequence data--fiducial, gain, and background files to load Olympus files, separate channels and save in matlab format
% Use this section ONLY if you have two-channel sequence data
% IMPORTANT - do not use this section for single snapshots, see below
for ii = 1:numel(allfiles)
    fprintf('%s:\n', allfiles{ii, 1});
    temp_vsi = bfopen(allfiles{ii,1});
    test = temp_vsi{1,1};
    filepath = cell2mat(test(1,2));
    filename = regexprep(filepath, ';.*$', '');
    [~, name, ext] = fileparts(filename);
    filename = {name};
    Channel_1 = [];
    Channel_2 = [];
%     sequence = [];
    for jj = 1:numel(test(:,1))
        if isEven(jj) == 1
            Channel_1 = cat(3,Channel_1,test{jj,1}); %even frames (488)
        else
            Channel_2 = cat(3,Channel_2,test{jj,1}); %odd frames  (561)
        end    
    end
%     sequence = Channel_1;
%     save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'sequence');
%     sequence = Channel_2;
%     save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'sequence');  

    save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
    save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
    % fig_temp1 = dipshow(Channel_1);
    % saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'tif')
    % saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'fig')
    % fig_temp2 = dipshow(Channel_2);
    % saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'tif')
    % saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'fig')
    close all
end
fprintf('Done 2-channel sequence data.\n');



%% Section B) Use this if you have 2-channel, single snapshot data--SiMPull 488/647
% IMPORTANT - Only use this section for single snapshots
for ii = 1:numel(allfiles)
    fprintf('%s:\n', allfiles{ii, 1});
    temp_vsi = bfopen(allfiles{ii,1});
    test = temp_vsi{1,1};
    filepath = cell2mat(test(1,2));
    filename = regexprep(filepath, ';.*$', '');
    [~, name, ext] = fileparts(filename);
    filename = {name};
    Channel_1 = [];
    Channel_2 = [];
    Channel_Merge = [];
    if length(test) == 1
        Channel_1 = test{1,1};
        save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
    elseif length(test) == 2
        Channel_1 = test{2,1};
        Channel_2 = test{1,1};
        save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
        save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
    elseif length(test) == 3
        Channel_1 = test{2,1};
        Channel_2 = test{1,1};
        Channel_Merge = test{3,1};
        save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
        save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
        % save(fullfile(SaveDir,[filename{1,1} '_Channel_Merge']), 'Channel_Merge');
    end
    if exist('Channel_1')
        fig_temp1 = dipshow(Channel_1);
        saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'tif')
        saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'fig')
    end
    if exist('Channel_2')
        fig_temp2 = dipshow(Channel_2);
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'tif')
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'fig')
    end

% Make merged image 
    red=Channel_2;
    green=Channel_1;
    co=joinchannels('RGB', stretch(red), stretch(green), stretch(red));
    fig_temp3 = dipshow(co);
    saveas(fig_temp3,fullfile(SaveDir,[filename{1,1} '_Merged']),'tif')
    saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Merged']),'fig')

    if exist('Channel_Merge')
        fig_temp2 = dipshow(Channel_2);
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_Merge']),'tif')
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_Merge']),'fig')
    end
    close all
end
fprintf('Done 1-channel snapshot data.\n');

%% Section C) 2-channel sequence data taken sequentially, ie 488/561
% SiMPull data where each channel has its own vsi fileSequential sequence.
%For two channel data taken sequentially, there is an empty frame. In the
%561 data, taken first so, the 1st image/channel
%should be kept and labeled as Channel 1. 
% Using the file list, it looks at if the file appears at an even or odd
% place on the list. I did this for consistency so it did not factor in the
% number on the file name. The odd places should always be 561, or Channel 1.
% By taking the even numbers and taking channel two from those, you get
% the 488/GFP data.  
% Use this section ONLY if you have two-channel sequence data
% IMPORTANT - Pay attention to channel numbers in the joining channels
% section. this might be opposite from the FR and GFP channel assignments. 

%this goes through the file list and looks at if it appears as even or odd in
%the sequence

for ii = 1:numel(allfiles)
    if rem(ii,2) == 0 %if there is no remainder after dividing ii by 2, is even
        fprintf('%s:\n', allfiles{ii, 1});
        temp_vsi = bfopen(allfiles{ii,1});
        test = temp_vsi{1,1};
        filepath = cell2mat(test(1,2));
        filename = regexprep(filepath, ';.*$', '');
        [~, name, ext] = fileparts(filename);
        filename = {name};
        Channel_2 = [];
        for jj = 1:numel(test(:,1))
            if isEven(jj) == 1 %if the frame is even, it is the GFP data
                Channel_2 = cat(3,Channel_2,test{jj,1}); 
            end
        end
        save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
        fig_temp1 = dipshow(Channel_2);
        saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'tif')
        saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'fig')
        close all
        
    elseif rem(ii,2) == 1 %if there is a remainder, it is odd
        fprintf('%s:\n', allfiles{ii, 1});
        temp_vsi = bfopen(allfiles{ii,1});
        test = temp_vsi{1,1};
        filepath = cell2mat(test(1,2));
        filename = regexprep(filepath, ';.*$', '');
        [~, name, ext] = fileparts(filename);
        filename = {name};
        Channel_1 = [];
        for jj = 1:numel(test(:,1))
            if isEven(jj)==0 %if the frame is odd it is the mCherry data
                Channel_1 = cat(3,Channel_1,test{jj,1});
            end 
        end
        save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
        fig_temp2 = dipshow(Channel_1);
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'tif')
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'fig')
        close all
   
%     sequence = Channel_1;
%     save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'sequence');
%     sequence = Channel_2;
%     save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'sequence');  

    % save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
    % save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
    % fig_temp1 = dipshow(Channel_1);
    % saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'tif')
    % saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'fig')
    %     fig_temp2 = dipshow(Channel_2);
    % saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'tif')
    % saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'fig')
    % close all
    end
end
fprintf('Done 2-channel sequence data taken sequentially.\n');
%% SECTION D) Use this if you have single channel sequence data
% IMPORTANT - do not use this section for single snapshots, see below
% This section has not been fully tested on single channel sequence data as of
% 06/15/2021
for ii = 1:numel(allfiles)
    fprintf('%s:\n', allfiles{ii, 1});
    temp_vsi = bfopen(allfiles{ii,1});
    test = temp_vsi{1,1};
    filepath = cell2mat(test(1,2));
    filename = regexprep(filepath, ';.*$', '');
    [~, name, ext] = fileparts(filename);
    filename = {name};
    Channel_1 = cat(3,test{:,1});
    sequence = Channel_1;
    save(fullfile(SaveDir,[filename{1,1} '_Channel_1']),'sequence');
end
fprintf('Done 1-channel sequence data.\n');