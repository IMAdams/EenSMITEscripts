% This script performs batch tracking from raw data, with the option to
% automatically apply (pre-prepared, external to this script) channel 
% registration.

%% Estimate the gain and offset (you can skip this and set manually in SMFs below)
% NOTE: This version only allows for a single calibration file for all
% data.

% Load the calibration files.
BeadFile = 'O:\Cell Path\Lidke Lab\Eric\Data\TIRF IX83\Emmy IX83 Calibration\Gain-512_20220311_1066_Channel_1.mat';
BgFile = 'O:\Cell Path\Lidke Lab\Eric\Data\TIRF IX83\Emmy IX83 Calibration\background-512_20220311_1070_Channel_1.mat';
load(BeadFile, 'Channel_1')
BeadData = single(Channel_1);
load(BgFile, 'Channel_1')
BgData = single(Channel_1);
clear Channel_1

% Use the DipImage function cal_readnoise() to estimate gain and offset.
CalResults = cal_readnoise(BeadData, BgData);
CameraGain = real(1 / CalResults(2)); % ADU/e-
CameraOffset = real(CalResults(4));
CameraReadNoise = real(CalResults(3));

%% Define directories, determine files to be tracked, and define parameters
% Define the path to the raw data.
RawDataDir = 'O:\Cell Path\Lidke Lab\IMAdams\Data\HeLa-ALFA-EGFR-KI\HelaALFAEGFR_10242025\data_dot_mat';
FilePattern = 'HeLa*.mat';

% Define a save sub-directory that will be placed under BaseDir.
ResultsDir = fullfile(RawDataDir, 'Results');

% Define the location of channel registration results.
TransformDir = 'empty';
TransformPattern = 'empty';

% Define the frame rate and pixel size.
FrameRate = 20; % fps
PixelSize = 0.067708; % micrometers

% Define our fitting/tracking parameters for each channel.
SMFChannel1 = smi_core.SingleMoleculeFitting;
SMFChannel1.Data.DataVariable = 'Channel_1';
SMFChannel1.Data.CameraType='SCMOS';
% SMFChannel1.Data.CameraGain = CameraGain; % ADU/e-
% SMFChannel1.Data.CameraOffset = CameraOffset;
% SMFChannel1.Data.CameraReadNoise = 9896.27;
SMFChannel1.Data.AnalysisID = 'Channel1';
SMFChannel1.Data.FrameRate = FrameRate;
SMFChannel1.Data.PixelSize = PixelSize;
SMFChannel1.Data.FileDir = RawDataDir;
SMFChannel1.Data.ResultsDir = ResultsDir;
SMFChannel1.BoxFinding.BoxSize = 11;
SMFChannel1.BoxFinding.BoxOverlap = 2;
SMFChannel1.BoxFinding.MinPhotons = 40;
SMFChannel1.Fitting.PSFSigma = 3;
SMFChannel1.Fitting.FitType = 'XYNBS';
SMFChannel1.Thresholding.AutoThreshLogL = false;
SMFChannel1.Thresholding.MinPValue = 0.000001;
SMFChannel1.Thresholding.MaxXY_SE = 0.75; % pixels
SMFChannel1.Thresholding.MinPhotons = 40;
SMFChannel1.Thresholding.MinPSFSigma = 0.75;
SMFChannel1.Thresholding.MaxPSFSigma = 4;
SMFChannel1.Tracking.D = 0.01; % pixel^2/frame
SMFChannel1.Tracking.K_on = 0.9; % 1/frame
SMFChannel1.Tracking.K_off = 0.1; % 1/frame
SMFChannel1.Tracking.MaxDistFF = 10; % pixels
SMFChannel1.Tracking.MaxDistGC = 10; % pixels
SMFChannel1.Tracking.MaxFrameGap = 10; % frames
SMFChannel1.Tracking.MinTrackLength = 5;
SMFChannel1.Tracking.TrajwiseD = true;
SMFChannel1.Tracking.NIterMaxBatch = 10;
SMFChannel1.Tracking.NIterMax = 1; % maybe safer to use 1 and change NIterMaxBatch
SMFChannel2 = copy(SMFChannel1);
SMFChannel2.Data.DataVariable = 'Channel_2';
SMFChannel2.Data.AnalysisID = 'Channel2';
% SMFChannel2.Data.DataROI = [1, 129, 128, 256, 1, 1];
SMFChannel2.Tracking.D = 0.01; % pixel^2/frame
SMFChannel2.Tracking.K_on = 0.9; % 1/frame
SMFChannel2.Tracking.K_off = 0.1; % 1/frame

% Prepare instances of the smi.SPT class for each channel and set some
% parameters.
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
SPTChannel2 = smi.SPT(SMFChannel2, 0);
SPTChannel2.GenerateMovies = false;
SPTChannel2.MovieParams.UnitFlag = 1;
SPTChannel2.MovieParams.Resolution = 0; % dpi, but 0 uses default screen resolution
SPTChannel2.GeneratePlots = true;
SPTChannel2.IsTestRun = false;
SPTChannel2.TransformDir = TransformDir;
SPTChannel2.TransformPattern = TransformPattern;
SPTChannel2.FindFiles = true;
SPTChannel2.FilePattern = FilePattern;
SPTChannel2.Verbose = 1;

%% Loop through the eligible raw data files and perform the tracking. 
[TRChannel1, SMDChannel1, SMDPreThreshCh1, FilesChannel1, TransformsChannel1] = ...
    SPTChannel1.batchTrack();
[TRChannel2, SMDChannel2, SMDPreThreshCh2, FilesChannel2, TransformsChannel2] = ...
    SPTChannel2.batchTrack();
fprintf('Batch tracking complete.\n')

%% Display a random movie from each channel.
% Define some movie parameters.
RandomFileIndex = randi(numel(FilesChannel1));
DisplayParams.AutoCrop = 0;
DisplayParams.UnitFlag = 0;
DisplayParams.MaxTrajLength = 40;

% Prepare the channel 1 movie.
LD = smi_core.LoadData;
[FilePathCh1, RandomFileCh1] = fileparts(FilesChannel1{RandomFileIndex});
SMFForMovieCh1 = copy(SMFChannel1);
SMFForMovieCh1.Data.FileDir = FilePathCh1;
SMFForMovieCh1.Data.FileName = [RandomFileCh1, '.mat'];
[~, RawDataChannel1] = LD.loadRawData(SMFForMovieCh1, 1, ...
    SMFForMovieCh1.Data.DataVariable);
load(fullfile(SPTChannel1.SMF.Data.ResultsDir, ...
    [RandomFileCh1, '_Channel1_Results.mat']), 'TR', 'SMD', 'SMDPreThresh')
MovieMaker = smi_vis.GenerateMovies;
MovieMaker.TR = TR;
MovieMaker.SMD = SMD;
% MovieMaker.SMD = SMDPreThresh; % allows us to look at thresholded locs.
MovieMaker.RawData = RawDataChannel1;
MovieMaker.SMF = SMFForMovieCh1;
MovieMaker.Params = DisplayParams;
MovieMaker.gui()
MovieMaker.GUIFigure.Name = 'Channel 1';

% Prepare the channel 2 movie.
LD = smi_core.LoadData;
[FilePathCh2, RandomFileCh2] = fileparts(FilesChannel2{RandomFileIndex});
SMFForMovieCh2 = copy(SMFChannel2);
SMFForMovieCh2.Data.FileDir = FilePathCh2;
SMFForMovieCh2.Data.FileName = [RandomFileCh2, '.mat'];
[~, RawDataChannel2] = LD.loadRawData(SMFForMovieCh2, 1, ...
    SMFForMovieCh2.Data.DataVariable);
load(fullfile(SPTChannel2.SMF.Data.ResultsDir, ...
    [RandomFileCh2, '_Channel2_Results.mat']), 'TR', 'SMD', 'SMDPreThresh')
MovieMaker = smi_vis.GenerateMovies;
MovieMaker.TR = TR;
MovieMaker.SMD = SMD;
% MovieMaker.SMD = SMDPreThresh; % allows us to look at thresholded locs.
MovieMaker.RawData = RawDataChannel2;
MovieMaker.SMF = SMFForMovieCh2;
MovieMaker.Params = DisplayParams;
MovieMaker.gui()
MovieMaker.GUIFigure.Name = ...
    'Channel 2 (channel reg. applied to tracks but not raw data)';
