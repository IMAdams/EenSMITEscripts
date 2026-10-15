
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
    fullfile(start_DataDir, 'Amphiregulin-HER2\20260825_RBL-EGFR-HER2\01_data_dot_mat')
    fullfile(start_DataDir, 'Amphiregulin-HER2\20260813_RBL_EGFR-HER2\01_data_dot_mat')
    };
filePatterns = {... % order with rawDataDirs
    '*RBL*Channel_1*.mat'...
    '*RBL*Channel_1*.mat'...
    % '*EGFR*Channel_2*.mat'...
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
% 
% % % Define our fitting/tracking parameters for each channel.
SMFChannel1 = smi_core.SingleMoleculeFitting();
% SMFChannel1.Data.CameraGain = gain; % ADU/e-
% SMFChannel1.Data.CameraOffset = offset;
% SMFChannel1.Data.CameraReadNoise = 0;

SMFChannel1.Data.AnalysisID = char(strcat(analysisID, "ch1"));
SMFChannel1.Data.DataVariable = 'Channel_1';
% SMFChannel1.Data.DataROI = [1,513,512,1024];
SMFChannel1.Data.FrameRate = FrameRate;
SMFChannel1.Data.PixelSize = PixelSize;
SMFChannel1.Data.FileDir = RawDataDir;
SMFChannel1.Data.ResultsDir = ResultsDir; 
SMFChannel1.Data.CameraType = 'SCMOS';
SMFChannel1.Data.CalibrationFilePath = 'C:\Users\imadams\Documents\smite workspace\iX83_gainoffset\CH1_GainCal_256.mat';
SMFChannel1.BoxFinding.BoxSize = 9;
SMFChannel1.BoxFinding.BoxOverlap = 5;
SMFChannel1.BoxFinding.MinPhotons = 85;
SMFChannel1.Fitting.PSFSigma = 1.8;
SMFChannel1.Fitting.FitType = 'XYNBS';
SMFChannel1.Thresholding.On = true;
SMFChannel1.Thresholding.MaxXY_SE = 0.8; % pixels
SMFChannel1.Thresholding.MinPValue = 0.01;
% SMFChannel1.Thresholding.AutoThreshLogL = true;
SMFChannel1.Thresholding.MinPSFSigma = 1.1;
SMFChannel1.Thresholding.MaxPSFSigma = 2.3;
SMFChannel1.Thresholding.MinPhotons = 50;
SMFChannel1.Thresholding.MaxBg = 8;
% SMFChannel1.FrameConnection.On = true;
SMFChannel1.Tracking.D = 0.01; % pixel^2/frame
SMFChannel1.Tracking.K_on = 0.9; % 1/frame
SMFChannel1.Tracking.K_off = 0.1; % 1/frame
SMFChannel1.Tracking.MaxDistFF = 5; % pixels
SMFChannel1.Tracking.MaxDistGC = 8; % pixels
SMFChannel1.Tracking.MaxFrameGap = 2; % frames
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
SPTChannel1.TransformDir =RawDataDir;
SPTChannel1.TransformPattern = TransformPattern;
SPTChannel1.FindFiles = true;
SPTChannel1.FilePattern = FilePattern;
SPTChannel1.Verbose = 1;
% 
[TRChannel1, SMDChannel1, SMDPreThreshCh1, FilesChannel1, TransformsChannel1] = ...
   SPTChannel1.batchTrack();
fprintf('Batch tracking complete.\n')

% this uses a local function I made to convert the SMF to a .txt. 
fprintf('Batch tracking complete.\n')
save(fullfile(ResultsDir,"SMFChannel1.mat"), "SMFChannel1")
filename = "SMFChannel1";
fileIn = fullfile(ResultsDir, strcat(filename, ".mat"));
fileOut = fullfile(ResultsDir, strcat(filename, ".txt"));
SMF = load(fileIn);
SMF = SMF.SMFChannel1;
IMA_exportStructToText(SMF, fileOut)

end