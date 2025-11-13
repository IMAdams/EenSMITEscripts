% This script performs batch tracking from raw data, with the option to
% automatically apply (pre-prepared, external to this script) channel 
% registration.

% loops over multiple days of analysis. 

rawDataDirs = {...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250404_cos7_alfaHER2'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250410_cos7ALFAHER2_antiALFA_EGFqd'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250417_Cos-7-ALFA-HER2'...
'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250418_Cos7-ALFA-HER2\Test'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250422_cos7ALFAHER2'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250425_cos7ALFAHER2'...
%'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20241212_cos7ALFAHER2_egfQD'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250402_cos7ALFAHER2lowSort2'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\240426-Cos7_ALFA-EGFR_EGFtwocolor'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\20240515_alfaEGFR_EGF_twocolor'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\20240516_alfaEGFR_EGF_twocolor'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\20240517_alfaEGFR_EGF_twocolor'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\20240621_alfaEGFR_EGF_twocolor'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\20240816_alfa_EGFR_EGF_twocolor'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\20240820_alfa_EGFR_EGF_twocolor'...
% 'O:\Cell Path\Lidke Lab\IMAdams\Data\Cos-7 ALFA-EGFR\SPT iX71 alfaEGFR cos7 cells\20240821_alfa_EGFR_EGF_2color'...
   }; 

 % filePattern in order with rawDataDirs
filePatterns = {...
% '*_cos*.mat'...
% '*_cos*.mat'...
% '*_Cos*.mat'...
'*_Cos*.mat'...
% '*_Cos*.mat'...
% '*_Cos*.mat'...
% '*cos7*.mat'...
% '*cos7*.mat'...
% '*_Cos*.mat'...
% '*_cos7*.mat'...
% '*_alfa*.mat'...
% '*_alfa*.mat'...
% '*_alfa*.mat'...
% '*_alfa*.mat'...
% '*_alfa*.mat'...
% '*_alfa*.mat'...
};

try
    % Assuming cellArray1 and cellArray2 are your cell arrays
    numElements1 = numel(rawDataDirs);
    numElements2 = numel(filePatterns);

    if numElements1 == numElements2
        fprintf('Both cell arrays have %d elements.\n', numElements1);
    else
        fprintf('Cell array 1 has %d elements, cell array 2 has %d elements.\n', ...
            numElements1, numElements2);
    end

catch ME
    fprintf('Error: %s\n', ME.message);
    fprintf('Error occurred in: %s\n', ME.identifier);
end

for i = 1:numel(rawDataDirs)

% Load the calibration files.
RawDataDir = rawDataDirs{i};
BeadsFileName = smi_helpers.getFileNames(RawDataDir, 'Bead*.mat');
BeadFile = fullfile(RawDataDir, BeadsFileName);
BgFileName = smi_helpers.getFileNames(RawDataDir, 'Background*.mat');
BgFile = fullfile(RawDataDir, BgFileName);

load(BeadFile{1}, 'sequence');
BeadData = single(sequence);
load(BgFile{1}, 'sequence');
BgData = single(sequence);
clear sequence

% Use the DipImage function cal_readnoise() to estimate gain and offset.
CalResults = cal_readnoise(BeadData, BgData);
close all
CameraGain = real(1 / CalResults(2)); % ADU/e-
CameraOffset = real(CalResults(4));
CameraReadNoise = CalResults(3);

RawDataDir = rawDataDirs{i};
FilePattern = filePatterns{i};

% Define a save sub-directory that will be placed under BaseDir.
ResultsDir = fullfile(RawDataDir, 'Results1020');

% Define the location of channel registration results.
TransformPattern = 'RegistrationTransform*.mat'; % this is the output of channelreg script
TransformDir = RawDataDir;

% Define the frame rate and pixel size.
FrameRate = 20; % fps
PixelSize = 0.1667; % micrometers
  
% Define our fitting/tracking parameters for each channel.
SMFChannel1 = smi_core.SingleMoleculeFitting();
% SMF.Data
SMFChannel1.Data.FileDir = RawDataDir;
SMFChannel1.Data.ResultsDir = ResultsDir;
SMFChannel1.Data.AnalysisID = 'Channel1';
SMFChannel1.Data.CameraGain = CameraGain; % ADU/e-
SMFChannel1.Data.CameraOffset = CameraOffset;
% Running with camera readnoise
SMFChannel1.Data.CameraReadNoise = CameraReadNoise;
SMFChannel1.Data.DataROI = [1, 1, 128, 128, 1, 1];
SMFChannel1.Data.FrameRate = FrameRate;
SMFChannel1.Data.PixelSize = PixelSize;
% SMF.BoxFinding
SMFChannel1.BoxFinding.BoxSize = 7;
SMFChannel1.BoxFinding.BoxOverlap = 3;
SMFChannel1.BoxFinding.MinPhotons = 100;
% SMF.Fitting
SMFChannel1.Fitting.PSFSigma = 1.1811;
SMFChannel1.Fitting.FitType = 'XYNBS';
% SMF.Thresholding
SMFChannel1.Thresholding.On = true;
SMFChannel1.Thresholding.MaxXY_SE = 0.25; % pixels
SMFChannel1.Thresholding.MinPValue = 0.0001;
SMFChannel1.Thresholding.AutoThreshLogL = false;
% SMFChannel1.Thresholding.AutoThreshPrctile = 1e-05;
SMFChannel1.Thresholding.MinPSFSigma = 0.6;
SMFChannel1.Thresholding.MaxPSFSigma = 2;
SMFChannel1.Thresholding.MinPhotons = 100;
% ensure frame connect is on
SMFChannel1.FrameConnection.On = true;
% SMF.Tracking
SMFChannel1.Tracking.D = 0.01; % pixel^2/frame
SMFChannel1.Tracking.TrajwiseD = true;
SMFChannel1.Tracking.K_on = 0.9; % 1/frame
SMFChannel1.Tracking.K_off = 0.1; % 1/frame
SMFChannel1.Tracking.MaxDistFF = 2; % pixels
SMFChannel1.Tracking.MaxDistGC = 5; % pixels
SMFChannel1.Tracking.MaxFrameGap = 50; % frames
SMFChannel1.Tracking.MinTrackLength = 10;
SMFChannel1.Tracking.NIterMaxBatch = 10;
SMFChannel1.Tracking.NIterMax = 1; % maybe safer to use 1 and change NIterMaxBatch
SMFChannel2 = copy(SMFChannel1);
SMFChannel2.Data.AnalysisID = 'Channel2';
SMFChannel2.Data.DataROI = [1, 129, 128, 256, 1, 1];
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
SPTChannel1.Verbose = 2;
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
SPTChannel2.Verbose = 2;

[TRChannel1, SMDChannel1, SMDPreThreshCh1, FilesChannel1, TransformsChannel1] = ...
   SPTChannel1.batchTrack();
% [TRChannel2, SMDChannel2, SMDPreThreshCh2, FilesChannel2, TransformsChannel2] = ...
%    SPTChannel2.batchTrack();
fprintf('Batch tracking complete.\n')

% save SMF for channels 1 and 2
save(fullfile(ResultsDir,"SMFChannel1.mat"), "SMFChannel1")
filename = "SMFChannel1";
fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
SMF = load(fileIn);
SMF = SMF.SMFChannel1;
% IMA_exportStructToText(SMF, fileOut)
% save(fullfile(ResultsDir,"SMFChannel2.mat"), "SMFChannel2")
% filename = "SMFChannel2";
% fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
% fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
% SMF = load(fileIn);
% SMF = SMF.SMFChannel2;
% IMA_exportStructToText(SMFChannel2, fileOut)
end
%% Display a random movie from each channel.
% Define some movie parameters.
% RandomFileIndex = randi(numel(FilesChannel1));
DisplayParams.AutoCrop = 0;
DisplayParams.UnitFlag = 0;
DisplayParams.MaxTrajLength = 40;
% 
% % Prepare the channel 1 movie.
% LD = smi_core.LoadData;
% [FilePathCh1, RandomFileCh1] = fileparts(FilesChannel1{RandomFileIndex});
% SMFForMovieCh1 = copy(SMFChannel1);
% SMFForMovieCh1.Data.FileDir = FilePathCh1;
% SMFForMovieCh1.Data.FileName = [RandomFileCh1, '.mat'];
% [~, RawDataChannel1] = LD.loadRawData(SMFForMovieCh1, 1, ...
%     SMFForMovieCh1.Data.DataVariable);
% load(fullfile(SPTChannel1.SMF.Data.ResultsDir, ...
%     [RandomFileCh1, '_Channel1_Results.mat']), 'TR', 'SMD', 'SMDPreThresh')
MovieMaker = smi_vis.GenerateMovies;
MovieMaker.TR = TR;
MovieMaker.SMD = SMD;
% MovieMaker.SMD = SMDPreThresh; % allows us to look at thresholded locs.
MovieMaker.RawData = sequence(128:255,:,:);
MovieMaker.SMF = SMF;
MovieMaker.Params = DisplayParams;
MovieMaker.gui()
MovieMaker.GUIFigure.Name = 'Channel 2';

% Prepare the channel 2 movie.
% LD = smi_core.LoadData;
% [FilePathCh2, RandomFileCh2] = fileparts(FilesChannel2{RandomFileIndex});
% SMFForMovieCh2 = copy(SMFChannel2);
% SMFForMovieCh2.Data.FileDir = FilePathCh2;
% SMFForMovieCh2.Data.FileName = [RandomFileCh2, '.mat'];
% [~, RawDataChannel2] = LD.loadRawData(SMFForMovieCh2, 1, ...
%     SMFForMovieCh2.Data.DataVariable);
% load(fullfile(SPTChannel2.SMF.Data.ResultsDir, ...
%     [RandomFileCh2, '_Channel2_Results.mat']), 'TR', 'SMD', 'SMDPreThresh')
% MovieMaker = smi_vis.GenerateMovies;
% MovieMaker.TR = TR;
% MovieMaker.SMD = SMD;
% MovieMaker.SMD = SMDPreThresh; % allows us to look at thresholded locs.
% MovieMaker.RawData = RawDataChannel2;
% MovieMaker.SMF = SMFForMovieCh2;
% MovieMaker.Params = DisplayParams;
% MovieMaker.gui()
% MovieMaker.GUIFigure.Name = ...
%     'Channel 2 (channel reg. applied to tracks but not raw data)';
