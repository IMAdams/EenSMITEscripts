% my script for running HMM analysis
%% Reload the tracked results and search for dimer candidates.
% NOTE: 'FileDir' is a directory containing both channel 1 and channel 2
%        results files.  findDimerCandidatesFromFiles() will search this
%        directory for files matching FilePatterns{1} and treat those files
%        as channel 1 results.  It will then do the same for
%        FilePatterns{2} (channel 2 results) and then attempt to match the
%        channel 1 and channel 2 files based on their file names (e.g.,
%        Data01_Channel1_Results.mat would be paired to
%        Data01_Channel2_Results.mat).


DataDir = ('O:\Cell Path\Lidke Lab\IMAdams\Cos-7-ALFA-Her2-S310F\SPT ix71 Cos-7 ALFA-Her2 cells\20250410_cos7ALFAHER2_antiALFA_EGFqd\Results1');
ResultsDir = DataDir ;% fullfile(DataDir, 'Results');
FilePatterns = {'*Channel1_Results.mat'; '*Channel2_Results.mat'};
DimerSeparation = 0.27; % Separation between fluorophores on a dimer (pixels)

MaxSeparation = 16; % Max. separation for dimer candidates (pre-processing) (pixels)

DomainSeparation = 5; % Typical domain size for the free, dimer, domain model (pixels)
ModelSpecifier = 'DF';
    % 'DDF': dimer, domain, free
    % 'DF': dimer, free
MinValidPoints = 1; %
Verbose = 2;
[TRArray, SMFArray, FileList] = smi_stat.HMM.findDimerCandidatesFromFiles(...
ResultsDir, FilePatterns, DimerSeparation, MaxSeparation, MinValidPoints, Verbose);

%Add diffusion coefficients to the TRArray.

% NOTE: You can also do this on a per-trajectory basis if needed.

DC = 0.15; % Diffusion Coefficient 
TRArray(1).DiffusionCoefficient = DC;
TRArray(2).DiffusionCoefficient = DC;
saveDirectory= fullfile(ResultsDir, 'HMM_MVP10_Dimer0.27_MaxSep16DF');

%% Prepare the HMM class and run the analysis.
% SMF = smi_core.SingleMoleculeFitting;
% SMF.Data.FrameRate = 20; % fps, specific to the loaded TRs above
% SMF.Data.PixelSize = 0.1667; % micrometers, specific to the loaded TRs above
HMM = smi_stat.HMM(TRArray, SMFArray);
HMM.MaxSeparation = MaxSeparation;  % Max. separation for dimer candidates (pre-processing) (pixels)
HMM.DimerSeparation = DimerSeparation ; % Separation between fluorophores on a dimer (pixels)
HMM.DomainSeparation = DomainSeparation; % Typical domain size for the free, dimer, domain model (pixels)
HMM.ModelSpecifier = ModelSpecifier;
HMM.SaveDir = saveDirectory;
HMM.GeneratePlots = [true; true]; % [basic plots; summary plots]
HMM.GenerateMovies = true;
HMM.UnitFlag = 0; % 0 will plot in units of pixels
HMM.RateParametersGuess = []; % set to [] for model 'ddf'
HMM.Verbose = 3;
%HMM.StateNames = {'Dimer';'Domain';'Free'};
HMM.MovieParams.PercentileFloor = 10;
performFullAnalysis(HMM)
%% movies???
MovieParams = smi_vis.GenerateMovies.prepDefaults();
DefaultParams = smi_vis.GenerateMovies.prepDefaults();
DefaultParams.AutoCrop = true;
DefaultParams.IndicateDimer = true;
DefaultParams.IndicateDimerCandidate = true;
DefaultParams.TrajColors = [0, 255, 0; 255, 0, 0];
MovieParams = smi_helpers.padStruct(MovieParams, DefaultParams);
SaveDir = 'C:\Users\imadams\OneDrive - University of New Mexico Health Sciences Center\Desktop\20250123_qd625-alfaHEr2\movies';
if (~exist('SaveDir', 'var') || isempty(SaveDir))
    SaveDir = pwd();
end
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end
MovieParams.CropToDimerCandidates = 1;
%% 

load('C:\Users\imadams\OneDrive - University of New Mexico Health Sciences Center\Desktop\20250123_qd625-alfaHEr2\Results\HMM_Results\Condition1\Dimers\SMFArrayDimer.mat');
load('C:\Users\imadams\OneDrive - University of New Mexico Health Sciences Center\Desktop\20250123_qd625-alfaHEr2\Results\HMM_Results\Condition1\Dimers\TRArrayDimer.mat');
% Loop through dimer pairs and prepare the movies.
% for nn = 1:size(TRArray, 1)
nn = 25;
    % Load the raw data for this dimer pair, if possible.
    FullFile1 = fullfile(SMFArray(nn, 1).Data.FileDir, ...
        SMFArray(nn, 1).Data.FileName{1});
    if isfile(FullFile1)
        LD = smi_core.LoadData;
        [~, RawData1] = LD.loadRawData(SMFArray(nn, 1), 1, ...
            SMFArray(nn, 1).Data.DataVariable);
    else
        RawData1 = ones(TRArray(nn, 1).YSize, ...
            TRArray(nn, 1).XSize);
    end
    FullFile2 = fullfile(SMFArray(nn, 2).Data.FileDir, ...
        SMFArray(nn, 2).Data.FileName{1});
    if isfile(FullFile2)
        LD = smi_core.LoadData;
        [~, RawData2] = LD.loadRawData(SMFArray(nn, 2), 1, ...
            SMFArray(nn, 2).Data.DataVariable);
    else
        RawData2 = ones(TRArray(nn, 2).YSize, ...
            TRArray(nn, 2).XSize);
    end

    % Load the registration transform, if available (we'll use this to
    % transform the raw data so that it aligns nicely for the movie).
    if isfile(SMFArray(nn, 2).Data.RegistrationFilePath)
        % We'll always assume channel 2 is the registered channel!
        load(SMFArray(nn, 2).Data.RegistrationFilePath, ...
            'RegistrationTransform')
        RawData2 = smi_core.ChannelRegistration.transformImages(...
            RegistrationTransform{2}, RawData2);
    end

    % Rescale the raw data and create the false color overlay.
    RawData1 = smi_vis.contrastStretch(RawData1, [], ...
        MovieParams.PercentileCeiling, ...
        MovieParams.PercentileFloor, ...
        MovieParams.MinScaleIntensity);
    RawData2 = smi_vis.contrastStretch(RawData2, [], ...
        MovieParams.PercentileCeiling, ...
        MovieParams.PercentileFloor, ...
        MovieParams.MinScaleIntensity);
    RawDataRGB = ones([size(RawData1, 1:2), 3, size(RawData1, 3)]);
    MovieParams.RawDataColors = [128 128 128; 255 96 208];
    for cc = 1:3
        RawDataRGB(:, :, cc, :) = ...
            RawData1*MovieParams.RawDataColors(1, cc) ...
            + RawData2*MovieParams.RawDataColors(2, cc);
    end

    % Prepare and save the movie.
    MovieMaker = smi_vis.GenerateMovies(MovieParams);
    MovieMaker.RawData = RawDataRGB;
    MovieMaker.TR = TRArrayDimer(nn, :);
    MovieSavePath = fullfile(SaveDir, ...
        sprintf('DimerPair%i_movie.mp4', nn));
    MovieMaker.saveMovie(MovieSavePath)


%% 

% Loop through dimer pairs and prepare the movies.
for nn = 1:size(TRArray, 1)
    % Load the raw data for this dimer pair, if possible.
    FullFile1 = fullfile(SMFArray(nn, 1).Data.FileDir, ...
        SMFArray(nn, 1).Data.FileName{1});
    if isfile(FullFile1)
        LD = smi_core.LoadData;
        [~, RawData1] = LD.loadRawData(SMFArray(nn, 1), 1, ...
            SMFArray(nn, 1).Data.DataVariable);
    else
        RawData1 = ones(TRArray(nn, 1).YSize, ...
            TRArray(nn, 1).XSize);
    end
    FullFile2 = fullfile(SMFArray(nn, 2).Data.FileDir, ...
        SMFArray(nn, 2).Data.FileName{1});
    if isfile(FullFile2)
        LD = smi_core.LoadData;
        [~, RawData2] = LD.loadRawData(SMFArray(nn, 2), 1, ...
            SMFArray(nn, 2).Data.DataVariable);
    else
        RawData2 = ones(TRArray(nn, 2).YSize, ...
            TRArray(nn, 2).XSize);
    end

    % Load the registration transform, if available (we'll use this to
    % transform the raw data so that it aligns nicely for the movie).
    if isfile(SMFArray(nn, 2).Data.RegistrationFilePath)
        % We'll always assume channel 2 is the registered channel!
        load(SMFArray(nn, 2).Data.RegistrationFilePath, ...
            'RegistrationTransform')
        RawData2 = smi_core.ChannelRegistration.transformImages(...
            RegistrationTransform{2}, RawData2);
    end

    % Rescale the raw data and create the false color overlay.
    RawData1 = smi_vis.contrastStretch(RawData1, [], ...
        MovieParams.PercentileCeiling, ...
        MovieParams.PercentileFloor, ...
        MovieParams.MinScaleIntensity);
    RawData2 = smi_vis.contrastStretch(RawData2, [], ...
        MovieParams.PercentileCeiling, ...
        MovieParams.PercentileFloor, ...
        MovieParams.MinScaleIntensity);
    RawDataRGB = ones([size(RawData1, 1:2), 3, size(RawData1, 3)]);
    for cc = 1:3
        RawDataRGB(:, :, cc, :) = ...
            RawData1*MovieParams.RawDataColors(1, cc) ...
            + RawData2*MovieParams.RawDataColors(2, cc);
    end

    % Prepare and save the movie.
    MovieMaker = smi_vis.GenerateMovies(MovieParams);
    MovieMaker.RawData = RawDataRGB;
    MovieMaker.TR = TRArray(nn, :);
    MovieSavePath = fullfile(SaveDir, ...
        sprintf('DimerPair%i_movie.mp4', nn));
    MovieMaker.saveMovie(MovieSavePath)
end

    % Prepare and save the movie.
    MovieMaker = smi_vis.GenerateMovies(MovieParams);
    MovieMaker.RawData = RawDataRGB;
    MovieMaker.TR = TRArrayDimer(1, :);
    MovieSavePath = fullfile(SaveDir, ...
        sprintf('DimerPair%i_movie.mp4', 1));
    MovieMaker.saveMovie(MovieSavePath)


%%

f = figure()
[PlotFigure] = visualizeRegistrationResults(RegistrationTransform)