%% This is meant to convert single channel data into the Quadview format needed in the smite/SiMPull analysis.
%After running the 'Olympus_file_separation_EB' run this.

%FOR MICHAEL-the last two sections are for the registraton 488/561 and then
%for the SiMPull data taken sequentially. You can skip to there
start_DataDir = 'E:\LidkeLab\iX83\20251113_HeLa_ALFAEGFR_HATGFBR2\data_dot_mat';
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


%%
%20 frame fiducial, Gain, background data, 512x512
%ImSZ=256;
ImSZ=512;
% ImSZ=250;
% sequence = newim(1024,1024,3);
datasetJoined = newim(2*ImSZ,2*ImSZ,50);
% datasetJoined = newim(1024,1024,20);
%d atasetJoined = newim(512,512,20);
% datasetJoined = newim(500,500,20);
% SaveDir = fullfile(RawDataDir, 'fiducial');
% if not(isfolder(SaveDir))
%     mkdir(SaveDir);
% end
for ii=1:length(fileName1)
    DATA_1 = load([newFiledir, fileName1{ii}]);
    DATA_2 = load([newFiledir, fileName2{ii}]);
    
datasetJoined(ImSZ : 2*ImSZ - 1, 0 : ImSZ - 1, 0:end) = DATA_1.Channel_1; % top right

datasetJoined(0 : ImSZ - 1,0 : ImSZ - 1, 0:end) = DATA_2.Channel_2; % top left

    [pathstr,name,ext] = fileparts(fileName1{ii});
%     save([newFiledir 'joined_' name ext],'datasetJoined'); %change so name ext=fileNamecell
      save([newFiledir 'joined_' name ext],'datasetJoined');

end

% 
% Indx{1}{1} = ImSZ : 2*ImSZ - 1; Indx{1}{2} = 0 : ImSZ - 1;      % 488 [TR]
% Indx{2}{1} = 0 : ImSZ - 1;      Indx{2}{2} = 0 : ImSZ - 1;      % 561 [TL]
% Indx{3}{1} = 0 : ImSZ - 1;      Indx{3}{2} = ImSZ : 2*ImSZ - 1; % 642 [BL]

%%
%1 frame data,
%USed 6/9/23

ImSZ=256;
% ImSZ=512;
% ImSZ=250;
% sequence = newim(1024,1024,3);
datasetJoined = newim(2*ImSZ,2*ImSZ,4);

for ii=1:length(fileName1)
    DATA_1 = load([newFiledir, fileName1{ii}]);
    Channel_1=cat(3,zeros(256,256,1),DATA_1.Channel_1);
    Channel_1=cat(3,zeros(256,256,2),DATA_1.Channel_1);
    Channel_1=cat(3,zeros(256,256,3),DATA_1.Channel_1);
%     
    % Channel_1=cat(3,zeros(512,512,1),DATA_1.Channel_1);
    % Channel_1=cat(3,zeros(512,512,2),DATA_1.Channel_1);
    % Channel_1=cat(3,zeros(512,512,3),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(250,250,1),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(250,250,2),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(250,250,3),DATA_1.Channel_1);

%     Somehow data in 6/14-17 had 2 frames each so ned to cat 2 frames and exclude in
%     fileinfo.m
    DATA_2 = load([newFiledir, fileName2{ii}]);
    Channel_2=cat(3,zeros(256,256,1),DATA_2.Channel_2);%1024x1024 im
    Channel_2=cat(3,zeros(256,256,2),DATA_2.Channel_2);%1024x1024 im
    Channel_2=cat(3,zeros(256,256,3),DATA_2.Channel_2);%1024x1024 im

    % Channel_2=cat(3,zeros(512,512,1),DATA_2.Channel_2);%1024x1024 im
    % Channel_2=cat(3,zeros(512,512,2),DATA_2.Channel_2);%1024x1024 im
    % Channel_2=cat(3,zeros(512,512,3),DATA_2.Channel_2);%1024x1024 im
%     Channel_2=cat(3,zeros(250,250,1),DATA_2.Channel_2);%1024x1024 im
%     Channel_2=cat(3,zeros(250,250,2),DATA_2.Channel_2);%1024x1024 im
%     Channel_2=cat(3,zeros(250,250,3),DATA_2.Channel_2);%1024x1024 im

%     DATA_2 = DATA_2.Channel_2(:,:,1);%two frames here, getting rid of one for now
%     sequence(  512:1023,   0:511) = DATA_1.Channel_1; % 488
%     datasetJoined(  512:1023,   0:511, 0:end) = Channel_1; % 488 USE THIS FOR IX83
% datasetJoined(  250:499,   0:249, 0:end) = Channel_1; % 488


    datasetJoined(  256:511,   0:255, 0:end) = Channel_1; % 488
    % datasetJoined(512:1023, 0:511, 0:end) = Channel_1;
    % 
    % datasetJoined(  0:511, 512:1023, 0:end) = Channel_2; % 647 USE THIS FOR IX83


%     datasetJoined(  0:255, 0:255, 0:end) = Channel_2; % 561 USE THIS FOR IX83
    datasetJoined(  0:255, 256:511, 0:end) = Channel_2; % 647

%     datasetJoined(Indx{3}{1}, Indx{3}{2}, :) = ...
%           dataset(Indx{3}{1}, Indx{3}{2}, :, 1); % 635
%     datasetJoined(Indx{2}{1}, Indx{2}{2}, :) = ...
%           dataset(Indx{2}{1}, Indx{2}{2}, :, 2); % 561
%     datasetJoined(Indx{1}{1}, Indx{1}{2}, :) = ...
%           dataset(Indx{1}{1}, Indx{1}{2}, :, 3); % 488
    [pathstr,name,ext] = fileparts(fileName1{ii});
    save([newFiledir 'joined_' name ext],'datasetJoined'); %change so name ext=fileNamecell
end

%%
%1 frame data, 512x512 
% ImSZ=512;
ImSZ=256;

% sequence = newim(1024,1024,3);
datasetJoined = newim(2*ImSZ,2*ImSZ,4);

for ii=1:length(fileName1)
    DATA_1 = load([newFiledir, fileName1{ii}]);
    Channel_1=cat(3,zeros(256,256,1),DATA_1.Channel_1);
    Channel_1=cat(3,zeros(256,256,2),DATA_1.Channel_1);
    Channel_1=cat(3,zeros(256,256,3),DATA_1.Channel_1);
%     
%     Channel_1=cat(3,zeros(512,512,1),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(512,512,2),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(512,512,3),DATA_1.Channel_1);
%     Somehow data in 6/14-17 had 2 frames each so ned to cat 2 frames and exclude in
%     fileinfo.m
    DATA_2 = load([newFiledir, fileName2{ii}]);
    Channel_2=cat(3,zeros(256,256,1),DATA_2.Channel_2);%1024x1024 im
    Channel_2=cat(3,zeros(256,256,2),DATA_2.Channel_2);%1024x1024 im
    Channel_2=cat(3,zeros(256,256,3),DATA_2.Channel_2);%1024x1024 im

%     Channel_2=cat(3,zeros(512,512,1),DATA_2.Channel_2);%1024x1024 im
%     Channel_2=cat(3,zeros(512,512,2),DATA_2.Channel_2);%1024x1024 im
%     Channel_2=cat(3,zeros(512,512,3),DATA_2.Channel_2);%1024x1024 im
%     DATA_2 = DATA_2.Channel_2(:,:,1);%two frames here, getting rid of one for now
%     sequence(  512:1023,   0:511) = DATA_1.Channel_1; % 488
%     datasetJoined(  512:1023,   0:511, 0:end) = Channel_1; % 488 USE THIS FOR IX83
    datasetJoined(  256:511,   0:255, 0:end) = Channel_1; % 488
    % datasetJoined(511:1023, 0:511,)=Channel_1;
    % sequence(  0:511, 512:1023, 0:end) = DATA_2.Channel_2; % 647
%     datasetJoined(  0:511, 512:1023, 0:end) = Channel_2; % 647 USE THIS FOR IX83
%     sequence(  0:511, 512:1023) = DATA_2; % 647
%     datasetJoined(  0:255, 0:255, 0:end) = Channel_2; % 561 USE THIS FOR IX83
    datasetJoined(  0:255, 256:511, 0:end) = Channel_2; % 647

%     datasetJoined(Indx{3}{1}, Indx{3}{2}, :) = ...
%           dataset(Indx{3}{1}, Indx{3}{2}, :, 1); % 635
%     datasetJoined(Indx{2}{1}, Indx{2}{2}, :) = ...
%           dataset(Indx{2}{1}, Indx{2}{2}, :, 2); % 561
%     datasetJoined(Indx{1}{1}, Indx{1}{2}, :) = ...
%           dataset(Indx{1}{1}, Indx{1}{2}, :, 3); % 488
    [pathstr,name,ext] = fileparts(fileName1{ii});
    save([newFiledir 'joined_' name ext],'datasetJoined'); %change so name ext=fileNamecell
end


%% FOR MICHAEL
% 488/561 channel 20 frame fiducial, Gain, background data, 256x256
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



    [pathstr,name,ext] = fileparts(fileName1{ii});
%     save([newFiledir 'joined_' name ext],'datasetJoined'); %change so name ext=fileNamecell
      save([newFiledir 'joined_' name ext],'datasetJoined');

end

%% FOR MICHAEL, Data in new acquisition
%This is for data taken sequentially and assigned the red channel and Ch1
%and green adn Ch 2. The location on the joined will still have the green
%in the top right and the red in top left)
%20 frame fiducial, Gain, background data, 512x512
ImSZ=256;

datasetJoined = newim(2*ImSZ,2*ImSZ,4);

for ii=1:length(fileName1)
    DATA_1 = load([newFiledir, fileName1{ii}]);
    DATA_1 = DATA_1.Channel_1(:,:,1);%two frames here, getting rid of the second one
    Channel_1=cat(3,zeros(256,256,1),DATA_1);
    Channel_1=cat(3,zeros(256,256,2),DATA_1);
    Channel_1=cat(3,zeros(256,256,3),DATA_1);
%     
    % Channel_1=cat(3,zeros(512,512,1),DATA_1.Channel_1);
    % Channel_1=cat(3,zeros(512,512,2),DATA_1.Channel_1);
    % Channel_1=cat(3,zeros(512,512,3),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(250,250,1),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(250,250,2),DATA_1.Channel_1);
%     Channel_1=cat(3,zeros(250,250,3),DATA_1.Channel_1);


    DATA_2 = load([newFiledir, fileName2{ii}]);
    DATA_2 = DATA_2.Channel_2(:,:,1);%two frames here, getting rid of one 
    Channel_2=cat(3,zeros(256,256,1),DATA_2);
    Channel_2=cat(3,zeros(256,256,2),DATA_2);
    Channel_2=cat(3,zeros(256,256,3),DATA_2);

%     DATA_2 = DATA_2.Channel_2(:,:,1);%two frames here, getting rid of one for now
%     sequence(  512:1023,   0:511) = DATA_1.Channel_1; % 488
%     datasetJoined(  512:1023,   0:511, 0:end) = Channel_1; % 488 USE THIS FOR IX83
% datasetJoined(  250:499,   0:249, 0:end) = Channel_1; % 488


    datasetJoined(  256:511,   0:255, 0:end) = Channel_2; % 488


  datasetJoined(  0:255, 0:255, 0:end) = Channel_1; % 561 USE THIS FOR IX83
    

%     datasetJoined(Indx{3}{1}, Indx{3}{2}, :) = ...
%           dataset(Indx{3}{1}, Indx{3}{2}, :, 1); % 635
%     datasetJoined(Indx{2}{1}, Indx{2}{2}, :) = ...
%           dataset(Indx{2}{1}, Indx{2}{2}, :, 2); % 561
%     datasetJoined(Indx{1}{1}, Indx{1}{2}, :) = ...
%           dataset(Indx{1}{1}, Indx{1}{2}, :, 3); % 488
    [pathstr,name,ext] = fileparts(fileName1{ii});
    save([newFiledir 'joined_' name ext],'datasetJoined'); %change so name ext=fileNamecell
end
