% This script generates channel registration transforms that can be used to
% transform SPT results.  When using this script, you want to keep
% ChannelReg.ManualCull = true and visually inspect the color coded plus
% and circle localizations overlain on the raw data.  Ideally, every
% emitter will be marked, and pairs of color coded pluses and circles
% should mark correct pairings between the two channels.  You can left
% click incorrect pairings to delete them from the transform.  Closing the
% figure will save the displayed localizations and compute the transform.

%% Create a list of the fiducial files.
DataDir = 'E:\LidkeLab\iX83\HelaALFAEGFR_10242025\VSI_to_mat';

FiducialFileDirs = {fullfile(DataDir)...
    %fullfile(DataDir, '1')...
  %fullfile(DataDir, '2')...
   %fullfile(DataDir, '3')...
 % fullfile(DataDir, '5')
 } ;
FiducialPattern = 'Fiduc*.mat';
% 
% % Create the Results directory within the obj.SaveDir (if needed).
% if (exist(obj.SMF.Data.ResultsDir, 'dir') ~= 7)
%     % exist() will return a 7 if the ResultsDir is an existing directory.
%     mkdir(obj.SMF.Data.ResultsDir)
% end
% Loop through all of the fiducial files and compute a transform.
SMF = smi_core.SingleMoleculeFitting;
SMF.Fitting.FitType = 'XYNBS';
SMF.Fitting.PSFSigma = 1.3;
SMF.BoxFinding.MinPhotons = 50;
% SMF.Thresholding.On = false; 
SMF.Thresholding.MinPValue = 0.01;
ChannelReg = smi_core.ChannelRegistration();
ChannelReg.ManualCull = true;
ChannelReg.SplitFormat = [1, 2];
ChannelReg.NNeighborPoints = 12;
ChannelReg.SeparationThreshold = 10;
ChannelReg.TransformationType = 'lwm';
ChannelReg.AutoscaleFiducials = true;
for ff = 1:numel(FiducialFileDirs)
    FiducialFileStruct = dir(fullfile(FiducialFileDirs{ff}, FiducialPattern));
    FiducialFileNames = {FiducialFileStruct.name}.';
    SMF.Data.FileDir = FiducialFileDirs{ff};
    ChannelReg.SMF = SMF;
    NFiducials = numel(FiducialFileNames);
    for ii = 1:NFiducials
        ChannelReg.SMF.Data.FileName = FiducialFileNames(ii);
        ChannelReg.findTransform();
        ChannelReg.exportTransform(FiducialFileDirs{ff});
        fprintf('Transform %i of %i computed: %s\n', ...
            ii, NFiducials, FiducialFileNames{ii})
    end
% end

% look at the transform and reg error if desired. 

% get RMSE, RMSE-LOO
FixedCoordinates = ChannelReg.Coordinates{1}(:,:,1);
MovingCoordinates = ChannelReg.Coordinates{2}(:,:,2);
RMSE = sqrt(mean(ChannelReg.estimateRegistrationError(...
    ChannelReg.RegistrationTransform{2}, ...
    MovingCoordinates, FixedCoordinates)));
% RMSELOO = sqrt(mean(ChannelReg.estimateRegErrorLOO(...
%     ChannelReg.TransformationType, {ChannelReg.NNeighborPoints}, ...
%     MovingCoordinates, FixedCoordinates)));

% Visualize the performance of the channel registration.
FixedImages = ChannelReg.FiducialImages(:, :, 1);
MovingImages = ChannelReg.FiducialImages(:, :, 2);
PlotFigure = figure();
ChannelReg.visualizeRegistrationResults(PlotFigure, ...
    ChannelReg.RegistrationTransform{2}, ...
    MovingCoordinates, FixedCoordinates, ...
    MovingImages, FixedImages);
saveas(gcf, fullfile(FiducialFileDirs{1}, 'CR1.png'));


% Visualize the registration error.
PlotFigure = figure();
PlotAxes = axes(PlotFigure);
ChannelReg.visualizeRegistrationError(PlotAxes, ...
    ChannelReg.RegistrationTransform{2}, ...
    MovingCoordinates, FixedCoordinates);
saveas(gcf, fullfile(FiducialFileDirs{1}, 'CR2.png'));

% Visualize the transform magnitude and gradient (this isn't usually useful
% unless something went very wrong, in which case it might be obvious in
% this plot!).
FiducialSize = [128, 128] + 1; % for the whole FOV 
%FiducialSize = diff(ChannelReg.FiducialROI([1, 2; 3, 4])) + 1;
PlotFigure = figure();
ChannelReg.visualizeCoordTransform(PlotFigure, ...
    ChannelReg.RegistrationTransform{2}, FiducialSize);
saveas(gcf, fullfile(FiducialFileDirs{1}, 'CR3.png'));


% Visualize the transforms effect on images.
% WARNING: This one hurts my eyes a bit... Also, it's not too useful unless
%          the transform is very dramatic.
PlotFigure = figure();
PlotAxes = axes(PlotFigure);
ChannelReg.visualizeImageTransform(PlotAxes, ...
    ChannelReg.RegistrationTransform{2}, FiducialSize);
saveas(gcf, fullfile(FiducialFileDirs{1}, 'CR4.png'));

end