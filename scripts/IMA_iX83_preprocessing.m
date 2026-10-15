%% Preprocessing two-channel data from iX83
% 1) converts from .vsi to .mat using% bio-formats 2) rejoins the separated
% channels into unified files. 3) compute the transforms. then I am running
% a powershell script to rename the files so they hav ethe time string in
% the filename, needed to match registration transforms to the data.


%% Directories and file set up
% I save the olympus data in it's own folder, then export the *.mat
% formatted into a new folder.

% base directory for processing
BaseDir = 'O:\Cell Path\Lidke Lab\IMAdams\Data\TGF-beta\SPT\20261008_MB231-ALFA-EGFR_2F9_tgfb1-QD655';

% Define the directory containing the raw data.
RawDataDir = ...
    '00_data_dot_vsi';
RawDataDir = fullfile(BaseDir, RawDataDir);

% Define a save directory in which the tracking results will be stored.

start_DataDir = '01_data_dot_mat';
start_DataDir = fullfile(BaseDir, start_DataDir);
if not(isfolder(start_DataDir))
    mkdir(start_DataDir);
end

allfiles = uipickfiles('FilterSpec', RawDataDir, 'REFilter', '\.vsi', ...
                       'Prompt', 'Select vsi Files to Combine')';


%%  1) converts from .vsi to .mat using% bio-formats 
for ii = 1:numel(allfiles) % !!
    fprintf('%s:\n', allfiles{ii, 1});
    temp_vsi = bfopen(allfiles{ii,1});
    test = temp_vsi{1,1};
    filepath = cell2mat(test(1,2));
    filename = regexprep(filepath, ';.*$', '');
    [~, name, ext] = fileparts(filename);
    filename = {name};
    Channel_1 = [];
    Channel_2 = [];
    for jj = 1:numel(test(:,1))
        if isEven(jj) == 1
            Channel_1 = cat(3,Channel_1,test{jj,1}); %even frames 
        else
            Channel_2 = cat(3,Channel_2,test{jj,1}); %odd frames  
        end    
    end

    save(fullfile(start_DataDir,[filename{1,1} '_Channel_1']), 'Channel_1');
    save(fullfile(start_DataDir,[filename{1,1} '_Channel_2']), 'Channel_2');
end

% for ii = 1:numel(allfiles)
%     fprintf('%s:\n', allfiles{ii, 1});
%     temp_vsi = bfopen(allfiles{ii,1});
%     test = temp_vsi{1,1};
%     filepath = cell2mat(test(1,2));
%     filename = regexprep(filepath, ';.*$', '');
%     [~, name, ext] = fileparts(filename);
%     filename = {name};
%     Channel_1 = [];
%     Channel_2 = [];
%     for jj = 1:numel(test(:,1))
%         if isEven(jj) == 1
%             Channel_1 = cat(3,Channel_1,test{jj,1}); %even frames 
%         else
%             Channel_2 = cat(3,Channel_2,test{jj,1}); %odd frames  
%         end    
%     end
% 
%     save(fullfile(start_DataDir,[filename{1,1} '_Channel_1']), 'Channel_1');
%     save(fullfile(start_DataDir,[filename{1,1} '_Channel_2']), 'Channel_2');
% end

%% Convert single channel data into the Quadview format
% TBD whether this is needed for all data files to use HMM. It is not
% needed for tracking non-transformed data. 

BaseDir = 'O:\Cell Path\Lidke Lab\IMAdams\Data\MDA-MB-231-ALFA-EGFR\SPT\20260707_MDA-ALFA-EGFR';
start_DataDir = '01_data_dot_mat';
start_DataDir = fullfile(BaseDir, start_DataDir);

% SiMPullUtils.joinSequentialChannels(start_DataDir, ImSZ, Indx);
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

joinFileDir = fullfile(BaseDir, '02_joined_dot_mat');
if not(isfolder(joinFileDir))
    mkdir(joinFileDir);
end

for ii=1:length(fileName1)
    DATA_1 = load([newFiledir, fileName1{ii}]);
    DATA_2 = load([newFiledir, fileName2{ii}]);
    
    ImSZ_y=length(DATA_1.Channel_1(:,1,1));
    ImSZ_x=length(DATA_1.Channel_1(1,:,1));
    ImSZ_z = length(DATA_1.Channel_1(1,1,:));

    datasetJoined = newim(2*ImSZ_x,2*ImSZ_y,ImSZ_z);

    datasetJoined(ImSZ_x : 2*ImSZ_x - 1, 0 : ImSZ_y - 1, 0:end) = DATA_1.Channel_1; % top right


    datasetJoined(0 : ImSZ_x - 1,0 : ImSZ_y - 1, 0:end) = DATA_2.Channel_2; % top left

    [pathstr,name,ext] = fileparts(fileName1{ii});
    save(fullfile(joinFileDir, ['joined_' name ext]),'datasetJoined', '-v7.3');

end

%% 
% Set to false to not use the results of channel registration.
UseRegistration = true;

% Set to true if using the Olympus IX83 (producing vsi files).
Olympus = true;

% Used for channel registration:
% Channel numbers that are being used in this analysis.  Note that Channel 1
% must always be present, and either 2 or 3 (or both) should also be present.
%Channels = [1, 2, 3];
%Channels = [1, 3];
Channels = [1, 2];
% User defined labels for the channels.
ChannelLabel = {'Qdot655', 'Qdot585'};
ImSZ_x = 512;
ImSZ_y = 512;
Indx{1}{1} = ImSZ_x : 2*ImSZ_x - 1; Indx{1}{2} = 0 : ImSZ_y - 1;      % [TR]
Indx{2}{1} = 0 : ImSZ_x - 1;      Indx{2}{2} = 0 : ImSZ_y - 1;      % [TL]
Indx{3}{1} = 0 : ImSZ_x - 1;      Indx{3}{2} = ImSZ_y : 2*ImSZ_y - 1; % [BL]

%% SimPull function to get gain and offset

[gain, offset] = SiMPullUtils.findGainOffset(joinFileDir);

%% ---------- Channel Registration ----------
% Perform registration between channels 2 through max(Channels) with respect to
% fixed channel 1.  The nanogrid files used for computing the geometric
% transforms between channels should be stored in DataDir/fiducial.  These
% consist of images of beads set on rectangular grids, to which a locally
% weighted mean transformation will be computed to align the moving grids
% (channels 2 or greater) to the fixed grid in channel 1.  This process will
% produce files GeomTransform_*.mat saved in in the same directory that
% contains the original fiducial data files.  For example, GeomTransfrom_FR_642
% contains the transformation that maps channel 3 onto channel 1.
% RegDataChannel*.mat files are temporary files that are created to properly
% call the SMITE regsitration routines.

fidDir='O:\Cell Path\Lidke Lab\IMAdams\Data\HeLa-ALFA-EGFR-KI\SPT\20260129_ALFAEGFR-HATGFBRII\02_joined_dot_mat';

fiducFiles = uipickfiles('FilterSpec', fidDir, 'REFilter', '\.mat', ...
                       'Prompt', 'Select fiducial Files to Combine')';
NNeighborPoints = 10;   % must be at least 6 but 6 sometimes gives errors

for i = 1:numel(fiducFiles)
  SiMPullUtils.channelRegistration(fiducFiles{i}, Indx, Channels, ...
                                    ChannelLabel, NNeighborPoints);
end
ImSZ_y=512
ImSZ_x=512
% Number of nearest neighbor points used to compute the geometric transforms
% between channels.
NNeighborPoints = 10;   % must be at least 6 but 6 sometimes gives errors
% if Olympus
%   SiMPullUtils.channelRegOlympus(start_DataDir, Indx, Channels, ...
% %                                  ChannelLabel, NNeighborPoints);
% % else
%    SiMPullUtils.channelRegistration(start_DataDir, Indx, Channels, ...
%                                     ChannelLabel, NNeighborPoints);
% % end


