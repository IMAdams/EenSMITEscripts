%% Define directories, determine files to be analyzed and separated
%%%IMPORTANT: I separated this script to process the registration, gain, background data from the Olympus IX83%%%
%%%When you choose files, only choose appropriate files in groups for the sections below

% Define the directory containing the raw data.

RawDataDir = ...
    'C:\Users\imadams\Documents\smite workspace\SPTtesttracking';

allfiles = uipickfiles('FilterSpec', RawDataDir, 'REFilter', '\.vsi', ...
                       'Prompt', 'Select vsi Files to Combine')';

% Define a save directory in which the tracking results will be stored.
SaveDir = fullfile(RawDataDir, 'Fiddotmat');
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
%    

%%  2 Channel Sequence data--fiducial, gain, and background files to load Olympus files, separate channels and save in matlab format
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
    fig_temp1 = dipshow(Channel_1);
    saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'tif')
    saveas(fig_temp1,fullfile(SaveDir,[filename{1,1} '_Channel_1']),'fig')
    fig_temp2 = dipshow(Channel_2);
    saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'tif')
    saveas(fig_temp2,fullfile(SaveDir,[filename{1,1} '_Channel_2']),'fig')
    close all
end
fprintf('Done 2-channel sequence data.\n');


%% Join data in QuadView--This is meant to convert single channel data into the Quadview format needed in the smite/SiMPull analysis.
%After running the 'Olympus_file_separation_EB'(above) run this.You will
%select all the channel 1 files and then the Channel 2

start_DataDir = SaveDir;
filedir = uigetdir(start_DataDir, 'Directory containing files to be tracked');
% savedir=fullfile(filedir,'QuadView');
% if ~exist(savedir, 'dir')
%   mkdir(savedir);
% end

[fileName1, newFiledir]=uigetfile(fullfile(start_DataDir, '*_1.mat'), ...
                                 'Select sequential channels',     ...
                                 'MultiSelect','on');
                             
[fileName2, newFiledir]=uigetfile(fullfile(start_DataDir, '*_2.mat'), ...
                                 'Select sequential channels',     ...
                                 'MultiSelect','on');
if ~iscell(fileName1)
   fileName1 = {fileName1};
end

if ~iscell(fileName2)
   fileName2 = {fileName2};
end

%Might need to change here depending on data length
% sequence = newim(1024,1024,20);
% sequence = newim(1024,1024,3); 
% sequence = newim(1024,1024);
% sequence = newim(512,512,20);
%% FOR MICHAEL
% 488/561 channel 20 frame fiducial, Gain, background data, 256x256
%To do 488/647, just uncomment the 647 line and comment out the 561 line. 
ImSZ=256;
% ImSZ=512;
% ImSZ=250;
% sequence = newim(1024,1024,3);
% datasetJoined = newim(2*ImSZ,2*ImSZ,4);
% datasetJoined = newim(1024,1024,20);
datasetJoined = newim(512,512,20);
% datasetJoined = newim(500,500,20);
% SaveDir = fullfile(RawDataDir, 'fiducial');
% if not(isfolder(SaveDir))
%     mkdir(SaveDir);
% end
for ii=1:length(fileName1)
    DATA_1 = load([newFiledir, fileName1{ii}]);
    DATA_2 = load([newFiledir, fileName2{ii}]);
    

    
    datasetJoined(  256:511,   0:255, 0:end) = DATA_1.Channel_1; % 488

    datasetJoined(   0:255, 0:255, 0:end) = DATA_2.Channel_2; % 561 USE THIS FOR IX83
   % datasetJoined(  0:255, 256:511, 0:end) = DATA_2.Channel_2; % 647 USE
   % THIS FOR IX83-might need to change to Channel 3 depending on how
   % Micheal changes



    [pathstr,name,ext] = fileparts(fileName1{ii});
%     save([newFiledir 'joined_' name ext],'datasetJoined'); %change so name ext=fileNamecell
      save([newFiledir 'joined_' name ext],'datasetJoined');

end