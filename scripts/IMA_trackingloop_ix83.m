% This script performs batch tracking from raw data collected on the
% TIRF scope (Olympus iX83). Files need to be separated using the file
% separation script. This is bootstrapped to use the EMCCD gain calibration
% code, rather than the sCMOS code, which has bug for SPT. 

% processed data converted from .vsi format to .mat.
baseDir = "E:\LidkeLab\iX83\20251113_HeLa_ALFAEGFR_HATGFBR2\data_dot_mat";

% [gain, offset] = SiMPullUtils.findGainOffset(baseDir);
BeadFile = fullfile(baseDir,'joined_Beads_20251113_12176_Channel_1.mat');
BgFile = fullfile(baseDir,'joined_Background_20251113_12179_Channel_1.mat');
load(BeadFile, "datasetJoined")
BeadData = single(datasetJoined);
load(BgFile, 'datasetJoined')
BgData = single(datasetJoined);
clear datasetJoined

% Use the DipImage function cal_readnoise() to estimate gain and offset.
CalResults = cal_readnoise(BeadData, BgData);
% close all;
CameraGain = real(1 / CalResults(2)); % ADU/e-
CameraOffset = real(CalResults(4));
CameraReadNoise = real(CalResults(3));
%%
% Define the frame rate and pixel size.
FrameRate = 20; % fps
PixelSize = 0.0541; % micrometers 
ImSZ = 512;
analysisID = '1114';

% Define the location of channel registration results.
% Haven't tried to make yet - Ian 
TransformPattern = 'GeomTransform*.mat';
TransformDir = baseDir;

rawDataDirs = {...
    fullfile(baseDir)...
        fullfile(baseDir)...
    };
filePatterns = {... % order with rawDataDirs
    'HeLa*Channel_1.mat'...
    'HeLa*Channel_2.mat'...
    };

resultsDir = strcat("Results_", analysisID);
%% track channel 1

RawDataDir = rawDataDirs{1};
FilePattern = filePatterns{1};

% Define a save sub-directory that will be placed under BaseDir.
ResultsDir = fullfile(RawDataDir, resultsDir);


% Define our fitting/tracking parameters for each channel.
SMFChannel1 = smi_core.SingleMoleculeFitting();
SMFChannel1.Data.CameraGain = CameraGain; % ADU/e-
SMFChannel1.Data.CameraOffset = CameraOffset;
SMFChannel1.Data.CameraReadNoise = CameraReadNoise;
SMFChannel1.Data.AnalysisID = char(strcat(analysisID, "ch1"));
SMFChannel1.Data.DataVariable = 'Channel_1';
SMFChannel1.Data.DataROI = [ImSZ : 2*ImSZ- 1, 0 : ImSZ - 1];
SMFChannel1.Data.FrameRate = FrameRate;
SMFChannel1.Data.PixelSize = PixelSize;
SMFChannel1.Data.FileDir = RawDataDir;
SMFChannel1.Data.ResultsDir = ResultsDir;
SMFChannel1.BoxFinding.BoxSize = 13;
SMFChannel1.BoxFinding.BoxOverlap = 5;
SMFChannel1.BoxFinding.MinPhotons = 25;
SMFChannel1.Fitting.PSFSigma = 1.8;
SMFChannel1.Fitting.FitType = 'XYNBS';
SMFChannel1.Thresholding.On = true;
SMFChannel1.Thresholding.MaxXY_SE = 0.35; % pixels
SMFChannel1.Thresholding.MinPValue = 0.19;
SMFChannel1.Thresholding.AutoThreshLogL = false;
SMFChannel1.Thresholding.MinPSFSigma = 1.0;
SMFChannel1.Thresholding.MaxPSFSigma = 3.0;
SMFChannel1.Thresholding.MinPhotons = 50;
SMFChannel1.FrameConnection.On = true;
SMFChannel1.Tracking.D = 0.01; % pixel^2/frame
SMFChannel1.Tracking.K_on = 0.9; % 1/frame
SMFChannel1.Tracking.K_off = 0.1; % 1/frame
SMFChannel1.Tracking.MaxDistFF = 5; % pixels
SMFChannel1.Tracking.MaxDistGC = 20; % pixels
SMFChannel1.Tracking.MaxFrameGap = 10; % frames
SMFChannel1.Tracking.MinTrackLength = 10;
SMFChannel1.Tracking.TrajwiseD = true;
SMFChannel1.Tracking.NIterMaxBatch = 10;
SMFChannel1.Tracking.NIterMax = 1; % maybe safer to use 1 and change NIterMaxBatch

% Prepare instance of the smi.SPT class 
SPTChannel1 = smi.SPT(SMFChannel1, 0);
SPTChannel1.GenerateMovies = false;
SPTChannel1.MovieParams.UnitFlag = 1;
SPTChannel1.MovieParams.Resolution = 0; % dpi, but 0 uses default screen resolution
SPTChannel1.GeneratePlots = true;
SPTChannel1.IsTestRun = false;
SPTChannel1.TransformDir = TransformDir;
SPTChannel1.TransformPattern = TransformPattern;
SPTChannel1.FindFiles = true;
SPTChannel1.FilePattern = FilePattern;
SPTChannel1.Verbose = 1;

[TRChannel1, SMDChannel1, SMDPreThreshCh1, FilesChannel1, TransformsChannel1] = ...
   SPTChannel1.batchTrack();
fprintf('Batch tracking complete.\n')

% % this uses a local function I made to convert the SMF to a .txt. 
% fprintf('Batch tracking complete.\n')
% save(fullfile(ResultsDir,"SMFChannel1.mat"), "SMFChannel1")
% filename = "SMFChannel1";
% fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
% fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
% SMF = load(fileIn);
% SMF = SMF.SMFChannel1;
% IMA_exportStructToText(SMF, fileOut)
%% track channel 2

RawDataDir = rawDataDirs{2};
FilePattern = filePatterns{2};

% Define a save sub-directory that will be placed under BaseDir.
ResultsDir = fullfile(RawDataDir, resultsDir);
ImSZ = 512;

% Define our fitting/tracking parameters for each channel.
SMFChannel2 = smi_core.SingleMoleculeFitting();
SMFChannel2.Data.CameraGain = CameraGain; % ADU/e-
SMFChannel2.Data.CameraOffset = CameraOffset;
SMFChannel2.Data.CameraReadNoise = CameraReadNoise;
SMFChannel2.Data.AnalysisID = char(strcat(analysisID, "ch2"));
SMFChannel2.Data.DataVariable = 'Channel_2';
% SMFChannel2.Data.DataROI = [0 : ImSZ - 1,0 : ImSZ - 1];
SMFChannel2.Data.FrameRate = FrameRate;
SMFChannel2.Data.PixelSize = PixelSize;
SMFChannel2.Data.FileDir = RawDataDir;
SMFChannel2.Data.ResultsDir = ResultsDir;
SMFChannel2.BoxFinding.BoxSize = 13;
SMFChannel2.BoxFinding.BoxOverlap = 5;
SMFChannel2.BoxFinding.MinPhotons = 25;
SMFChannel2.Fitting.PSFSigma = 1.8;
SMFChannel2.Fitting.FitType = 'XYNBS';
SMFChannel2.Thresholding.On = true;
SMFChannel2.Thresholding.MaxXY_SE = 0.35; % pixels
SMFChannel2.Thresholding.MinPValue = 0.19;
SMFChannel2.Thresholding.AutoThreshLogL = false;
SMFChannel2.Thresholding.MinPSFSigma = 1.0;
SMFChannel2.Thresholding.MaxPSFSigma = 3.0;
SMFChannel2.Thresholding.MinPhotons = 50;
SMFChannel2.FrameConnection.On = true;
SMFChannel2.Tracking.D = 0.01; % pixel^2/frame
SMFChannel2.Tracking.K_on = 0.9; % 1/frame
SMFChannel2.Tracking.K_off = 0.1; % 1/frame
SMFChannel2.Tracking.MaxDistFF = 5; % pixels
SMFChannel2.Tracking.MaxDistGC = 20; % pixels
SMFChannel2.Tracking.MaxFrameGap = 10; % frames
SMFChannel2.Tracking.MinTrackLength = 10;
SMFChannel2.Tracking.TrajwiseD = true;
SMFChannel2.Tracking.NIterMaxBatch = 10;
SMFChannel2.Tracking.NIterMax = 1; % maybe safer to use 1 and change NIterMaxBatch

% Prepare instance of the smi.SPT class 
SPTChannel2 = smi.SPT(SMFChannel2, 0);
SPTChannel2.GenerateMovies = false;
SPTChannel2.MovieParams.UnitFlag = 1;
SPTChannel2.MovieParams.Resolution = 0; % dpi, but 0 uses default screen resolution
SPTChannel2.GeneratePlots = true;
SPTChannel2.IsTestRun = false;
% SPTChannel2.TransformDir = TransformDir;
% SPTChannel2.TransformPattern = TransformPattern;
SPTChannel2.FindFiles = true;
SPTChannel2.FilePattern = FilePattern;
SPTChannel2.Verbose = 1;

[TRChannel2, SMDChannel2, SMDPreThreshCh2, FilesChannel2, TransformsChannel2] = ...
   SPTChannel2.batchTrack();
fprintf('Channel 2 tracking complete.\n')

% save(fullfile(ResultsDir,"SMFChannel2.mat"), "SMFChannel2")
% filename = "SMFChannel2";
% fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
% fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
% SMF = load(fileIn);
% SMF = SMF.SMFChannel2;
% IMA_exportStructToText(SMF, fileOut)
