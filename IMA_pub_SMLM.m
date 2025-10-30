% IMA_publishSMLM script to test publish scripting

%% Define the list of CoverslipDirs to analyze, create the SMF object and load an optimized SMF

% SCMOS Calibration filepath
SCMOSDir = 'C:\Users\imadams\Documents\smite workspace\SMLM_optimize\GainCalibration-2024_04_11_16_54_01.mat';

% Define base directory for analysis
baseDir = 'C:\Users\imadams\Documents\smite workspace\20250822_sorted_CHOHAEGFR\Export_data\Export_data_test_fits';

% Define coverslips directories relative to baseDir
coverslipDirs = {...
    'L858R_EGF\SMF_Loop'...
    'WT Rest\SMF_Loop'...
    };

% Define SMF files relative to baseDir
SMF_Files = {
    'SMF_DIAG_LOOP\SMF_1020_A.mat'...
    'SMF_DIAG_LOOP\SMF_1020_B.mat'...
    'SMF_DIAG_LOOP\SMF_1020_C.mat'...
    };

% SMF variables to be encoded in save directory name
keyVars = {...
    'BoxFinding.BoxSize'...
    'BoxFinding.BoxOverlap'...
    };

% Function to get value from nested struct using string path
function val = getNestedField(s, pathStr)
fields = strsplit(pathStr, '.');
val = s;
for i = 1:length(fields)
    val = val.(fields{i});
end
end

% Define output directory relative to baseDir
resultsDir = 'SMF_DIAG_LOOP';

% create fullpaths and struct
coverslipDirs = fullfile(baseDir, coverslipDirs);
SMF_Files = fullfile(baseDir, SMF_Files);
SMF_Batch = struct('name', {}, 'path', {});
for i = 1:length(SMF_Files)
    [~, filename, ~] = fileparts(SMF_Files{i});
    SMF_Batch(i).name = filename;
    SMF_Batch(i).path = SMF_Files{i};
end

% main loop for iterating over SMF versions
for i = 1:length(SMF_Batch)
    fprintf('Loading %s from %s\n', SMF_Batch(i).name, SMF_Batch(i).path);
    load(SMF_Batch(i).path);
    SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMF);
    SMF.Data.CalibrationFilePath = SCMOSDir;
    SMF.Data.SEAdjust = 3 / (1000 * SMF.Data.PixelSize);

    % Loop through each CoverslipDir
    for ii = 1:length(coverslipDirs)
        currentDir = coverslipDirs{ii};
        Publish = smi.Publish(SMF);
        % Publish.CoverslipDir = CoverslipDir;
        Publish.CoverslipDir = currentDir;
        if ~isempty(keyVars)
            % Build parameter string with just the values
            paramStr = '';
            for j = 1:length(keyVars)
                fields = strsplit(keyVars{j}, '.');
                val = SMF;
                for i = 1:length(fields)
                val = val.(fields{i});
                end
                paramStr = [paramStr, '_', num2str(val)];
            end
            Publish.SaveBaseDir = fullfile(baseDir, resultsDir, ...
                ['Results_', SMF_Batch(1).name, paramStr]);
        else
            Publish.SaveBaseDir = fullfile(baseDir, resultsDir, ...
                ['Results_', SMF_Batch(1).name]);
        end
        Publish.CellList = [];
        Publish.Verbose = 1;
        Publish.GenerateSR = 1;
        Publish.GenerateImagingStats = 0; %
        Publish.GenerateOverlayStats = 0;
        Publish.ShiftToReg = 0;
        Publish.SRImageZoom = 20;
        sprintf('Processing directory %d of %d: %s\n', ii, length(coverslipDirs), currentDir)
        Publish.performFullAnalysis()
    end

end
