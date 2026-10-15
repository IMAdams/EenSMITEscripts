
%%iX83_trackingTwoChannels

% Aquisition parameters
% Define the frame rate and pixel size.
FrameRate = 20; % fps
PixelSize = 0.0677; % micrometers 
ImSZ = 256; % pixels length of image frame
% PixelSize = 0.0541; % micrometers 
% ImSZ = 512; % pixels length of image frame
analysisID = 'SMF082726';

% Directory containing joined*.mat preprocessed data
start_DataDir = 'O:\Cell Path\Lidke Lab\IMAdams\Data';

% Define the location of channel registration results.
TransformPattern = 'RegistrationTransform*.mat';
% TransformDir = start_DataDir;

rawDataDirs = {...
    fullfile(start_DataDir, 'Amphiregulin-HER2\20260813_RBL_EGFR-HER2\01_data_dot_mat')...
    % fullfile(start_DataDir, 'Amphiregulin-HER2\20260825_RBL-EGFR-HER2\01_data_dot_mat')...
    };
filePatterns = {... % order with rawDataDirs
    '*RBL*Channel_2*.mat'...
    % '*RBL*Channel_2*.mat'...
    % '*EGFR*Channel_2*.mat'...
    };
% calculate gain and offset
% [gain, offset] = SiMPullUtils.findGainOffset(start_DataDir);

resultsDir = strcat("Results_", analysisID);
set(0, 'DefaultFigureVisible', 'off');
dbstop if error;
%% track channel 2
 
for i = 1:numel(rawDataDirs)

RawDataDir = rawDataDirs{i};
FilePattern = filePatterns{i};

% Define a save sub-directory that will be placed under BaseDir.
ResultsDir = fullfile(RawDataDir, resultsDir);
if not(isfolder(ResultsDir))
    mkdir(ResultsDir);
end

% track channel 2

% % % Define our fitting/tracking parameters for each channel.
SMFChannel2 = smi_core.SingleMoleculeFitting();
% SMFChannel2.Data.CameraGain = gain; % ADU/e-
% SMFChannel2.Data.CameraOffset = offset;
% SMFChannel2.Data.CameraReadNoise = CameraReadNoise;
SMFChannel2.Data.AnalysisID = char(strcat(analysisID, "ch2"));
SMFChannel2.Data.DataVariable = 'Channel_2';
% SMFChannel2.Data.DataROI = [1,1,512,512];
SMFChannel2.Data.FrameRate = FrameRate;
SMFChannel2.Data.PixelSize = PixelSize;
SMFChannel2.Data.FileDir = RawDataDir;
SMFChannel2.Data.ResultsDir = ResultsDir;
SMFChannel2.Data.CameraType = 'SCMOS';
SMFChannel2.Data.CalibrationFilePath = 'C:\Users\imadams\Documents\smite workspace\iX83_gainoffset\CH2_GainCal_256.mat';
SMFChannel2.BoxFinding.BoxSize = 9;
SMFChannel2.BoxFinding.BoxOverlap = 5;
SMFChannel2.BoxFinding.MinPhotons = 85;
SMFChannel2.Fitting.PSFSigma = 1.8;
SMFChannel2.Fitting.FitType = 'XYNBS';
SMFChannel2.Thresholding.On = true;
SMFChannel2.Thresholding.MaxXY_SE = 0.8; % pixels
SMFChannel2.Thresholding.MinPValue = 0.01;
% SMFChannel2.Thresholding.AutoThreshLogL = true;
SMFChannel2.Thresholding.MinPSFSigma = 1.1;
SMFChannel2.Thresholding.MaxPSFSigma = 2.3;
SMFChannel2.Thresholding.MinPhotons = 25;
SMFChannel2.Thresholding.MaxBg = 8;
% SMFChannel2.Thresholding.InMeanMultiplier = inf; 
% SMFChannel2.FrameConnection.On = true;
SMFChannel2.Tracking.D = 0.01; % pixel^2/frame
SMFChannel2.Tracking.K_on = 0.9; % 1/frame
SMFChannel2.Tracking.K_off = 0.1; % 1/frame
SMFChannel2.Tracking.MaxDistFF = 5; % pixels
SMFChannel2.Tracking.MaxDistGC = 10; % pixels
SMFChannel2.Tracking.MaxFrameGap = 2; % frames
SMFChannel2.Tracking.MinTrackLength = 5;
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
SPTChannel2.TransformDir = RawDataDir;
SPTChannel2.TransformPattern = TransformPattern;
SPTChannel2.FindFiles = true;
SPTChannel2.FilePattern = FilePattern;
SPTChannel2.Verbose = 2;

[TRChannel2, SMDChannel2, SMDPreThreshCh2, FilesChannel2, TransformsChannel2] = ...
   SPTChannel2.batchTrack();
fprintf('Channel 2 tracking complete.\n')





save(fullfile(ResultsDir,"SMFChannel2.mat"), "SMFChannel2")
filename = "SMFChannel2";
fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
SMF = load(fileIn);
SMF = SMF.SMFChannel2;
IMA_exportStructToText(SMF, fileOut) 

end