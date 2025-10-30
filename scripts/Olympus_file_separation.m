%% Define directories, determine files to be analyzed and separated
%%%IMPORTANT: This script only works on 1 or 2 channel data from the Olympus IX83%%%
%%%When you choose files, only choose appropriate files in groups for the sections below

% Define the directory containing the raw data.
RawDataDir = 'E:\LidkeLab\iX83\HelaALFAEGFR_10242025';
allfiles = selectMultipleFiles(RawDataDir, '.vsi','Select Files to Combine');
% filterBy = '.vsi';
% allfiles = uipickfiles('FilterSpec', RawDataDir, ...
%         'REFilter', filterBy);
% Define a save directory in which the tracking results will be stored.
SaveDir = fullfile(RawDataDir, 'VSI_to_mat');
if not(isfolder(SaveDir))
    mkdir(SaveDir);
end
% if ~exist(fullfile(SaveDir))
%     mkdir(Savedir);
% end


%% Load Olympus files, separate channels and save in matlab format 
% Use this section ONLY if you have two-channel sequence data
% IMPORTANT - do not use this section for single snapshots, see below
for ii = 1:numel(allfiles)
    temp_vsi = bfopen(allfiles{ii,1});
    test = temp_vsi{1,1};
    filepath = cell2mat(test(1,2));
    filestr = strsplit(filepath,'\'); 
    filename_temp = cell2mat(filestr(end));
    filename = strsplit(filename_temp,'.vsi');
    Channel_1 = [];
    Channel_2 = [];
%     sequence = [];
%     for jj = 1:numel(test(:,1))
%          if isEven(jj) == 1
%             Channel_1 = cat(3,test{jj,1},Channel_1); %even frames
%         else
%             Channel_2 = cat(3,test{jj,1},Channel_2); %odd frames   
%          end    
%     end
    for jj = 1:numel(test(:,1))
         if isEven(jj) == 1
            Channel_1 = cat(3,Channel_1,test{jj,1}); %even frames
        else
            Channel_2 = cat(3,Channel_2,test{jj,1}); %odd frames   
         end    
    end
%     sequence = Channel_1;
%     save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'sequence');
%     sequence = Channel_2;
%     save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'sequence');    
    save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
    save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
    fig_temp1 = dipshow(Channel_1);
    saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'tif')
    saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'fig')
        fig_temp2 = dipshow(Channel_2);
    saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'tif')
    saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'fig')
    close all
end

%% Use this if you have single channel sequence data
% IMPORTANT - do not use this section for single snapshots, see below
% This section has not been fully tested on single channel sequence data as of 06/15/2021
for ii = 1:numel(allfiles)
    temp_vsi = bfopen(allfiles{ii,1});
    test = temp_vsi{1,1};
    filepath = cell2mat(test(1,2));
    filestr = strsplit(filepath,'\'); 
    filename_temp = cell2mat(filestr(end));
    filename = strsplit(filename_temp,'.vsi');
    Channel_1 = cat(3,test{:,1});
    sequence = Channel_1;
    save(fullfile(SaveDir,[filename{1,1} '_Channel_1']),'sequence');
end
%% Use this if you have single snapshot data
% IMPORTANT - Only use this section for single snapshots
for ii = 1:numel(allfiles)
    temp_vsi = bfopen(allfiles{ii,1});
    test = temp_vsi{1,1};
    filepath = cell2mat(test(1,2));
    filestr = strsplit(filepath,'\'); 
    filename_temp = cell2mat(filestr(end));
    filename = strsplit(filename_temp,'.vsi');
    Channel_1 = [];
    Channel_2 = [];
    Channel_Merge = [];
    if length(test) == 1
        Channel_1 = test{1,1};
        save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
%     elseif length(test) == 2
%         Channel_1 = test{1,1};
%         Channel_2 = test{2,1};
%         save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
%         save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
    else length(test) == 3
        Channel_1 = test{2,1};
        Channel_2 = test{1,1};
        Channel_Merge = test{3,1};
        save(fullfile(SaveDir,[filename{1,1} '_Channel_1']), 'Channel_1');
        save(fullfile(SaveDir,[filename{1,1} '_Channel_2']), 'Channel_2');
        save(fullfile(SaveDir,[filename{1,1} '_Channel_Merge']), 'Channel_Merge');
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
    if exist('Channel_Merge')
        fig_temp2 = dipshow(Channel_2);
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_Merge']),'tif')
        saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_Merge']),'fig')
    end
    close all
end
