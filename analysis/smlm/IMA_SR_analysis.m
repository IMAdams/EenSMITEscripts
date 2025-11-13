CoverslipDir = 'C:\Users\imadams\OneDrive - University of New Mexico Health Sciences Center\Desktop\DNA-PAINToptim\test_cells\Cell_03';
SCMOSDir = 'C:\Users\imadams\OneDrive - University of New Mexico Health Sciences Center\Desktop\SMLM_optimize\GainCalibration-2024_04_11_16_54_01.mat';
SMF = smi_core.SingleMoleculeFitting;
% load('/mnt/nas/cellpath/Genmab/Data/2024-04-23_HeLa_SaturatingIgG10min/20240501-400box6_MEadj_SMF_FC.mat');
SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMF);

SMF.Data.CalibrationFilePath = SCMOSDir;
SMF.Data.SEAdjust = 3 / (1000 * SMF.Data.PixelSize); % nm converted to pixels
%% adapted from example script from LidkeLab/smite/MATLAB/examples/Example_SMLM_script.m

Saving = true;

% --- 2D ---

% SMF = smi_core.SingleMoleculeFitting();

SMF.Data.ResultsDir = 'C:\Users\imadams\OneDrive - University of New Mexico Health Sciences Center\Desktop\DNA-PAINToptim\L858R_rest\Results';

if Saving
   if ~exist(SMF.Data.ResultsDir, 'dir')
      mkdir(SMF.Data.ResultsDir);
   end
end

SMF.Data.FileDir           = CoverslipDir;
   % File name (cell array of char array)
SMF.Data.FileName         = 'Data_2025-1-30-16-25-46.h5';
   % ID tagged onto saved results (char array)(Default='')
SMF.Data.AnalysisID       = 'een3';
   % 'EMCCD','SCMOS' (Default='EMCCD')
SMF.Data.CameraType       = 'SCMOS';
% set pixel size
SMF.Data.PixelSize = 0.097; %microns
SMF.Fitting.FitType = 'XYNBS';
SMF.Thresholding.MinPSFSigma = 1.1;
SMF.Thresholding.MaxPSFSigma = 1.5;

% Camera Gain, scalar or image (Default=1)
  
% SMF.Data.CameraGain       = 2.2;
   % Camera Offset, scalar or image (Default=0)
% SMF.Data.CameraOffset     = 100;
   % Perform thresholding? (Default=true)
SMF.Thresholding.On       = true;
   % Maximum allowed precision in x,y (Pixels)(Default=.2)
SMF.Thresholding.MaxXY_SE = 0.2;
   % Perform frame connection? (Default=true)
SMF.FrameConnection.On    = true;
   % Perform drift correction? (Default=true)
SMF.DriftCorrection.On    = true;
    % threshold based on the P value. (default=0.0100)
SMF.Thresholding.MinPValue = 0.0001;
% Nearest neighbor threshold (default = 0)
SMF.Thresholding.MinNumNeighbors = 2;
% playing with different sigma
SMF.Fitting.PSFSigma = 1.3; %default = 1
% playing with differnet MinPhotons
SMF.BoxFinding.MinPhotons = 100; % default is 200
% Create an SMLM object using the values in the SMF structure.
SMLMobj = smi.SMLM(SMF);
SMLMobj.Verbose = 5;
% Do a test fit, displaying all the results to the screen (if VerboseTest >= 5,
% otherwise saving the results in ResultsDir/TestFit).  Note that calling
% testFit from the smi.gui is equivalent to what we are doing here.
SMLMobj.VerboseTest = 5;
SMLMobj.testFit(1);
% Do a full analysis, saving results in ResultsDir.
SMLMobj.fullAnalysis();

% save the SMF
filename = fullfile(SMF.Data.ResultsDir,strcat(SMF.Data.AnalysisID,'SMF4.mat'));
save(filename, "SMF",'-mat')