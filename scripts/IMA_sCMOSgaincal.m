%% Create a gain calibration file for sCMOS camera
% Developed with Olympus iX83 with Hammamatsu Fusion sCMOS with Gemini
% W-View beamsplitter

% manually collected

Params = [];
MeanLevel = [];
VarLevel = [];
start_DataDir = 'O:\Cell Path\Lidke Lab\IMAdams\Data\HeLa-ALFA-EGFR-KI\SPT\20260107_SPT_HeLa_ALFAEGFR_HATGFBR2\GainCal010726\data_dot_mat';
[gainCalFiles, newFiledir]=uigetfile(fullfile(start_DataDir, '*joined*.mat'), ...
                                 'Select gain sequences',     ...
                                 'MultiSelect','on');
%%
for i = 1:numel(gainCalFiles)

    load([fullfile(newFiledir, gainCalFiles{i})], 'datasetJoined');
    gainSequence = single(datasetJoined);
    gainSequence = gainSequence (1:512, 1:512,:);
    MeanLevel = cat(3, MeanLevel, mean(gainSequence, 3));
    VarLevel = cat(3, VarLevel, var(single(gainSequence), [], 3));
 end
Params.MeanLevel = single(MeanLevel);
Params.VarLevel = single(VarLevel);

% Perform Least-Squares Fit on collected data:
% 
%     leastSquaresFit from smi_stat class performs a least squares (ls)
%     fit on the provided data.
%     This function computes a least squares fit for the paired data in
%     XData, YData, and Weights. Check
%     [`smi_stat.leastSquaresFit()`](../../MATLAB/+smi_stat/leastSquaresFit.m)
%     for further documentation.

    % Here [:, :, 1] is ls offset, [:, :, 2] is ls slope (Gain)
    Beta = NaN(size(MeanLevel, 1), size(MeanLevel, 2), 2);
    for ii = 1:size(MeanLevel, 1)
       disp(ii)
       for jj = 1:size(MeanLevel, 2)
          Beta(ii, jj, 1:2) = smi_stat.leastSquaresFit(                   ...
             squeeze(MeanLevel(ii, jj, :)), squeeze(VarLevel(ii, jj, :)), ...
             1 ./ squeeze(VarLevel(ii, jj, :)));
       end
    end

    dipshow(Beta(:, :, 2))
    dipshow(Beta(:, :, 1))
    figure;
    histogram(Beta(:, :, 2))
    xlabel('Gain (ADU/e-)')

    figure;
    histogram(MeanLevel(:,:,1))
    xlabel('ADU')
    Params.CCDVar=single(VarLevel(:,:,1));
    Params.Gain=single(Beta(:,:,2));
    Params.CCDOffset=single(MeanLevel(:,:,1));
    Params.CameraObj.ROI = [1, 1, 512, 512];
%Save the calibration:

    SaveDir = 'C:\Users\imadams\Documents\smite workspace\iX83_gainoffset';
    FileName = fullfile(SaveDir, ...
       ['GainCalibration-', smi_helpers.genTimeString(),'CH_2']);
    % Params.CameraObj.ROI = CameraSCMOS.ROI;
    % Params.LampPowerRange = single(LampPowerRange);
    save(FileName, 'Params', '-v7.3')

