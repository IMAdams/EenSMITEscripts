
%%iX83_trackingTwoChannels

% Aquisition parameters
% Define the frame rate and pixel size.
FrameRate = 20; % fps
PixelSize = 0.0541; % micrometers 
ImSZ = 512; % pixels length of image frame
analysisID = '031326';

% Directory containing joined*.mat preprocessed data
start_DataDir = 'O:\Cell Path\Lidke Lab\Angela\Data\IX83 QD-SPT\';

% Define the location of channel registration results.
TransformPattern = 'RegistrationTransform*.mat';
% TransformDir = start_DataDir;

rawDataDirs = {...
    fullfile(start_DataDir, '20250108_HEK-HA-RONmNG_sort2bMID_HIGH_Total\data_dot_mat')...
    };
filePatterns = {... % order with rawDataDirs
     '*RON*2.mat'...
    };
% calculate gain and offset
% [gain, offset] = SiMPullUtils.findGainOffset(start_DataDir);

resultsDir = strcat("Results_", analysisID);
set(0, 'DefaultFigureVisible', 'off');
dbstop if error;
%% track channel 1
 
for i = 1:numel(rawDataDirs)

RawDataDir = rawDataDirs{i};
FilePattern = filePatterns{i};

% Define a save sub-directory that will be placed under BaseDir.
ResultsDir = fullfile(RawDataDir, resultsDir);
if not(isfolder(ResultsDir))
    mkdir(ResultsDir);
end

% % Define our fitting/tracking parameters for each channel.
% SMFChannel1 = smi_core.SingleMoleculeFitting();
% % SMFChannel1.Data.CameraGain = gain; % ADU/e-
% % SMFChannel1.Data.CameraOffset = offset;
% % SMFChannel1.Data.CameraReadNoise = 0;
% 
% SMFChannel1.Data.AnalysisID = char(strcat(analysisID, "ch1"));
% SMFChannel1.Data.DataVariable = 'datasetJoined';
% SMFChannel1.Data.DataROI = []; % [1,513,512,1024];
% SMFChannel1.Data.FrameRate = FrameRate;
% SMFChannel1.Data.PixelSize = PixelSize;
% SMFChannel1.Data.FileDir = RawDataDir;
% SMFChannel1.Data.ResultsDir = ResultsDir; 
% SMFChannel1.Data.CameraType = 'SCMOS';
% SMFChannel1.Data.CalibrationFilePath = 'C:\Users\imadams\Documents\smite workspace\iX83_gainoffset\GainCalibration-2026_01_19_13_28_42CH_1.mat';
% SMFChannel1.BoxFinding.BoxSize = 11;
% SMFChannel1.BoxFinding.BoxOverlap = 3;
% SMFChannel1.BoxFinding.MinPhotons = 125;
% SMFChannel1.Fitting.PSFSigma = 1.8;
% SMFChannel1.Fitting.FitType = 'XYNBS';
% SMFChannel1.Thresholding.On = true;
% SMFChannel1.Thresholding.MaxXY_SE = 1.2; % pixels
% % SMFChannel1.Thresholding.MinPValue = 0.01;
% SMFChannel1.Thresholding.AutoThreshLogL = true;
% SMFChannel1.Thresholding.MinPSFSigma = 0.6;
% SMFChannel1.Thresholding.MaxPSFSigma = 4.0;
% SMFChannel1.Thresholding.MinPhotons = 100;
% SMFChannel1.FrameConnection.On = true;
% SMFChannel1.Tracking.D = 0.01; % pixel^2/frame
% SMFChannel1.Tracking.K_on = 0.9; % 1/frame
% SMFChannel1.Tracking.K_off = 0.1; % 1/frame
% SMFChannel1.Tracking.MaxDistFF = 5; % pixels
% SMFChannel1.Tracking.MaxDistGC = 10; % pixels
% SMFChannel1.Tracking.MaxFrameGap = 10; % frames
% SMFChannel1.Tracking.MinTrackLength = 10;
% SMFChannel1.Tracking.TrajwiseD = true;
% SMFChannel1.Tracking.NIterMaxBatch = 10;
% SMFChannel1.Tracking.NIterMax = 1; % maybe safer to use 1 and change NIterMaxBatch
% 
% % Prepare instance of the smi.SPT class 
% SPTChannel1 = smi.SPT(SMFChannel1, 0);
% SPTChannel1.GenerateMovies = false;
% SPTChannel1.MovieParams.UnitFlag = 1;
% SPTChannel1.MovieParams.Resolution = 0; % dpi, but 0 uses default screen resolution
% SPTChannel1.GeneratePlots = true;
% SPTChannel1.IsTestRun = false;
% SPTChannel1.TransformDir =RawDataDir;
% SPTChannel1.TransformPattern = TransformPattern;
% SPTChannel1.FindFiles = true;
% SPTChannel1.FilePattern = FilePattern;
% SPTChannel1.Verbose = 2;
% % 
% [TRChannel1, SMDChannel1, SMDPreThreshCh1, FilesChannel1, TransformsChannel1] = ...
%    SPTChannel1.batchTrack();
% fprintf('Batch tracking complete.\n')
% 
% % % this uses a local function I made to convert the SMF to a .txt. 
% fprintf('Batch tracking complete.\n')
% save(fullfile(ResultsDir,"SMFChannel1.mat"), "SMFChannel1")
% filename = "SMFChannel1";
% fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
% fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
% SMF = load(fileIn);
% SMF = SMF.SMFChannel1;
% IMA_exportStructToText(SMF, fileOut)

% track channel 2

% Define our fitting/tracking parameters for each channel.
SMFChannel2 = smi_core.SingleMoleculeFitting();
% SMFChannel2.Data.CameraGain = gain; % ADU/e-
% SMFChannel2.Data.CameraOffset = offset;
% SMFChannel2.Data.CameraReadNoise = CameraReadNoise;
SMFChannel2.Data.AnalysisID = char(strcat(analysisID, "ch2"));
SMFChannel2.Data.DataVariable = 'Channel_2';
SMFChannel2.Data.DataROI = [1,1,512,512];
SMFChannel2.Data.FrameRate = FrameRate;
SMFChannel2.Data.PixelSize = PixelSize;
SMFChannel2.Data.FileDir = RawDataDir;
SMFChannel2.Data.ResultsDir = ResultsDir;
SMFChannel2.Data.CameraType = 'SCMOS';
SMFChannel2.Data.CalibrationFilePath = 'C:\Users\imadams\Documents\smite workspace\iX83_gainoffset\GainCalibration-2026_01_19_14_04_43CH_2.mat';
SMFChannel2.BoxFinding.BoxSize = 11;
SMFChannel2.BoxFinding.BoxOverlap = 3;
SMFChannel2.BoxFinding.MinPhotons = 50;
SMFChannel2.Fitting.PSFSigma = 1.8;
SMFChannel2.Fitting.FitType = 'XYNBS';
SMFChannel2.Thresholding.On = true;
SMFChannel2.Thresholding.MaxXY_SE = 3.0; % pixels
% SMFChannel2.Thresholding.MinPValue = 0.01;
SMFChannel2.Thresholding.AutoThreshLogL = true;
SMFChannel2.Thresholding.MinPSFSigma = 1.0;
SMFChannel2.Thresholding.MaxPSFSigma = 3.0;
SMFChannel2.Thresholding.MinPhotons = 1;
SMFChannel2.Thresholding.MaxBg = inf;
SMFChannel2.Thresholding.InMeanMultiplier = 2.5; % trying 3/13/26 to filter out super bright stuff
SMFChannel2.FrameConnection.On = true;
SMFChannel2.Tracking.D = 0.01; % pixel^2/frame
SMFChannel2.Tracking.K_on = 0.9; % 1/frame
SMFChannel2.Tracking.K_off = 0.1; % 1/frame
SMFChannel2.Tracking.MaxDistFF = 5; % pixels
SMFChannel2.Tracking.MaxDistGC = 10; % pixels
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
SPTChannel2.TransformDir = RawDataDir;
SPTChannel2.TransformPattern = TransformPattern;
SPTChannel2.FindFiles = true;
SPTChannel2.FilePattern = FilePattern;
SPTChannel2.Verbose = 2;

[TRChannel2, SMDChannel2, SMDPreThreshCh2, FilesChannel2, TransformsChannel2] = ...
   SPTChannel2.batchTrack();
fprintf('Channel 2 tracking complete.\n')



Define our fitting/tracking parameters for each channel.

save(fullfile(ResultsDir,"SMFChannel2.mat"), "SMFChannel2")
filename = "SMFChannel2";
fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
SMF = load(fileIn);
SMF = SMF.SMFChannel2;
IMA_exportStructToText(SMF, fileOut) 

end