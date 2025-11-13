% Description:
% This toolbox concentrates on organizing SiMPull data, and then performing
% single molecule fitting measurements using SMITE.  Single molecule emitters
% are fit and represented by Gaussian blobs.  Emitters from separate channels
% are then considered to be colocalized if they are determined to be within one
% pixel of one another.  These overlapping events are then used in the final
% correction step to determine the percentage of overlap.  In this demo data,
% EGR-GFP is channel 1.  GFP fits are cross-referenced with the anti-pan
% pY-AF647 (channel 3) to determine the percentage of phosphorylation.  

% System requirements:
%
% MATLAB Toolboxes: Image Processing, Optimization, Parallel Computing,
%                   Statistics and Machine Learning
% DIPimage (https://diplib.org)
% SMITE
%
% SMITE version; @keith - how do we supply the appropriate SMITE version??

% Directories and files used/created by this toolbox:
%
% File structure:
%    DataDir            contains images to be tracked; location of joined
%                       files and results.mat/results.xls outputs
%       results.mat        contains contents of valuecell which tabulates fit
%                          counts and overlap percentages between all channels
%                          (see overlayFitting.m for details)
%        results.xls       Excel file containing this data for easy viewing
%    DataDir/fiducial   input fiducial,Gain, Background, and output GeomTransform* files
%       GainOffset*.png    plots produced by DIPimage cal_readnoise
%       GeomTransform*     geometric transforms produced by channel reg.
%       RegDataChannel*    temporary files needed for SMITE channel reg.
%    DataDir/totrack    contains various temporary info saved by overlayFitting
%       *_*.mat            raw data sequence files identified by channel
%       *_shiftedvals.mat  shifted coordinates of the given moving channel
%       *__tripleOverlapCoord.mat  triple overlap coordinates for 3 channels
%    DataDir/Figures    contains overlay images and summary plots
%       *_boxes.png        fit boxes (pre-thresholding) found for each channel
%       *_overlayed.png    overlay of the Gaussian blob image on the original
%                          sequence data before applying any shifts
%       *_overlayedchannels.*  overlayed Gaussian blob images of the 3 channels
%       *_overlayW*only.*  only show Gaussian blobs in the specified channel
%                          that overlay channel 1 blobs
%       *_raw.*            raw data presented in quad view (fig and png)
%       *_rawOverlayLoc.*  quad view of overlayed blobs wrt ch. 1 per channel
%       Hist_*_minPhotons.*   histogram of photons per localization for the
%                          specified value of Thresholding_MinPhotons in title
%    DataDir/Corrected  contains background corrected plots and tables
%       graph_*            percentage overlap bar graph with SEM error bars
%                          saved as both .fig and .png
%       resultsStats_*.xls Excel spreadsheet of mean fit values with
%                          corrections that are used in the bar graph
%       table_*.mat        temporary file used to create the Excel spreadsheet

%% ---------- Initialization ----------
% This section sets the data directory that will be used throughout the
% analysis.  It also sets the parameters used in the SiMPull calculations.  In
% summary, ImSZ defines         the image size, with a default of 256x256 pixels.  Indx
% defines the quadrant indexing scheme for data assumed compactified into a
% quadrant format (see below for details).  Channels specifies what channels
% are to be analyzed for colocalization, and ChannelLabels are names to be used
% for each channel.
%
% The demo data consists of datasets from three ROIs, each defining a separate
% image series per channel of data.  In this example, we care about channels 1
% and 3, containing 488 and 642 data, respectively.  Therefore, we set
% Channels = [1, 3]; to exclude the 561 data in channel 2.
%
% IMPORTANT NOTE: The Initialization section should always be run first after
% any new invocation of MATLAB to properly set things up.  The later sections
% depend on the values below, so it is important to properly define them
% whenever MATLAB is started afresh and SiMPull is being used.

% start_DataDir = 'O:\Cell Path\Lidke Lab\Emmy\Manuscripts\SiMPullMethods_Jove\SampleData_Test-ThreshPhot-475';
start_DataDir = '.';


ImSZ = 256;   % image dimension (pixel); assumed square
% Channel indexing in a (2 ImSZ) x (2 ImSZ) quad.  The data is organized into
% quadrants of a larger image.  Indx{i}{j} are the indices defining the ith
% channel (1, 2, 3) of the jth dimension (1 = x, 2 = y), so, for example, using
% the default image size (ImSZ) of 256, the indices are defined as follows:
%    Indx{1}{1} = 256 : 511;  Indx{1}{2} =   0 : 255;   % 488 [top right]
%    Indx{2}{1} =   0 : 255;  Indx{2}{2} =   0 : 255;   % 561 [top left]
%    Indx{3}{1} =   0 : 255;  Indx{3}{2} = 256 : 511;   % 642 [bottom left]
% Note that we are using DIPimage, so image coordinates are zero based.  To
% reference the top right quad corresponding to channel 1,
% quadImage(Indx{1}{1}, Indx{1}{2}); can be used.  The original data actually
% has a 3rd dimension indexed by frame number, so this indexing scheme is used
% for this data as sequence(Indx{1}{1}, Indx{1}{2}, lo_frame : hi_frame);.
% Moreover, for the fiducial images (see below), the frames are interleaved, so
% that the odd frame numbers (1:2:end) are one dataset while the even frames
% (2:2:end) are a second, requiring additional complexity in the indexing.

Indx{1}{1} = ImSZ : 2*ImSZ - 1; Indx{1}{2} = 0 : ImSZ - 1;      % 488 [TR]
Indx{2}{1} = 0 : ImSZ - 1;      Indx{2}{2} = 0 : ImSZ - 1;      % 561 [TL]
Indx{3}{1} = 0 : ImSZ - 1;      Indx{3}{2} = ImSZ : 2*ImSZ - 1; % 642 [BL]
% The definitions above are the only place where image size and quad selection
% indices are defined to simplify extending this code to other cameras.

% Specify movie files to be processed and optionally exclude bad frames via a
% fileinfo.m file located in the directory containing the files to be tracked
% in the format:
%
% dat(1).filename = 'nameWithoutExtension';
% % Define bad frames (start counting with 0 as DIPimage); leave empty or omit
% % entirely for all files if no bad frames present in the data.
% badFrames{1} = [];

% Set to false to not use the results of channel registration.
UseRegistration = true;

% Set to true if using the Olympus IX83 (producing vsi files).
Olympus = true;

% Used for channel registration:
% Channel numbers that are being used in this analysis.  Note that Channel 1
% must always be present, and either 2 or 3 (or both) should also be present.
%Channels = [1, 2, 3];
%Channels = [1, 3];
Channels = [1, 2];
% User defined labels for the channels.
ChannelLabel = {'G_488', 'R_561', 'FR_642'};

fprintf('Done initialization.\n');

%% ---------- Find Camera Gain and Offset ----------
%9-25-25-
% The demo Gain and Background files for this section are found in
% DataDir/fiducial. 
% Running this section will request an input image file (typically, a nanogrid
% beads image taken out of focus [see also Channel Registration below]) and a
% background/camera bias image.  The DIPimage function cal_readnoise will then
% be called.  A gain (in units of e-/ADU) and offset will be given in the first
% plot produced.  However, SMITE requires gains in ADU/e-, so the inverted
% value for gain as well as the value for the offset will be written to the
% console (and also returned to the variables 'gain' and 'offset' below).  The
% plots produced by DIPImage (GainOffset*.png) will be saved in
% DataDir/fiducials, where the input and background images should also
% reside.Upon running this section, select the Gain image sequence.
% Another file selection will pop up and direct you to select the
% Background image.  The gain and offset calculated from this demo data are
% used in the single molecule fitting section as:
% smf.Data_CameraGain = [17.9, 17.9, 17.9];
% smf.Data_CameraOffset = [112.6, 112.6, 112.6];
[gain, offset] = SiMPullUtils.findGainOffset(start_DataDir);

%% ---------- Channel Registration ----------
% Perform registration between channels 2 through max(Channels) with respect to
% fixed channel 1.  The nanogrid files used for computing the geometric
% transforms between channels should be stored in DataDir/fiducial.  These
% consist of images of beads set on rectangular grids, to which a locally
% weighted mean transformation will be computed to align the moving grids
% (channels 2 or greater) to the fixed grid in channel 1.  This process will
% produce files GeomTransform_*.mat saved in in the same directory that
% contains the original fiducial data files.  For example, GeomTransfrom_FR_642
% contains the transformation that maps channel 3 onto channel 1.
% RegDataChannel*.mat files are temporary files that are created to properly
% call the SMITE regsitration routines.

% Number of nearest neighbor points used to compute the geometric transforms
% between channels.
NNeighborPoints = 10;   % must be at least 6 but 6 sometimes gives errors
%if Olympus
%   SiMPullUtils.channelRegOlympus(start_DataDir, Indx, Channels, ...
%                                  ChannelLabel, NNeighborPoints);
%else
   SiMPullUtils.channelRegistration(start_DataDir, Indx, Channels, ...
                                    ChannelLabel, NNeighborPoints);
%end

%% ---------- Join Sequential Channels into a Quad Image ----------
%9-25-25- Suggestion-Ideally we should be able to use this for both SiMPull data and
%registration fiducials. I would like to put this as the second section so it
% can be ran top to bottom. However, this incorporates the empty frames in the
%quadview to work with the existing code. It would be ideal if we could get
%rid of the "padding" and just allow the data to have whatever z-dimension
%it has when read in. Additionally, this is set up for 488/561 data to
%occupy the two channels. It would be great if we could have the code look
%for either 561 or 647. This might require having the files be save in the
%first place to either be channel 2 for 561 or channel 3 for 647. If we do
%not care about having three channel analysis, we could scrap the quadview
%entirely and just have the channels be side by side to flip through.  

% This section consolidates the data from multiple channels into a quad view.
% The output is an image stack of all ROIs, presented in a quad view format
% with the green 488 channel in the top right quadrant and the far red 642
% channel in the bottom left.  The file is saved in the same directory as the
% original data, but with 'joined_' prepended to the filename.  Input is all
% datasets to be analyzed which the user :R!will select upon running this section.
% This combined file allows easy visualization for curating the ROIs in the
% next section and will only need to be run once. To view individual
% datasets, load the file and use dipshow(datasetJoined).  Use 'n' and
% 'p' to flip through the ROI. The first frame will be empty.

if Olympus
   SiMPullUtils.joinSeqChannelsOlympus(start_DataDir, ImSZ, Indx);
else
   SiMPullUtils.joinSequentialChannels(start_DataDir, ImSZ, Indx);
end

%% ---------- Remove Bad Frames and generate filename list ----------
% This section allows for visual curation of the data.  The raw data from the
% previous section has been organized to allow for easy review so that bad (out
% of focus or high density) images can be removed before analysis.  The goal of
% this section is to take the reorganized data from the previous section (files
% with the name 'joined_*'), review all the datasets at once, let the user
% determine the bad frames that should be excluded, and provide file naming
% outputs for use in 'fileinfo.m', which specifies what data should be fit in
% the next section.  A 4D DIP image holding all the datasets will be created,
% each having an extra empty frame.  After running, a figure window will
% appear.  Use 'n' and 'p' to flip through the frames for each dataset.  Use
% 'f' and 'b' to go between the datasets. Note that the first dataset consists 
% entirely of empty frames.  Assess which frames should be excluded.  
% On the command console, file names will be printed that are
% formatted to list the frames that need to be excluded in a separate variable
% called badFrames.  Copy the file names and fill in the badFrames{...} in
% 'fileinfo.m'.  The provided fileinfo.m is formatted for the example data set.
% In this demo data, we exclude frame 3 frome dataset 1,frames 2 and 3 from
% dataset 4, all frames from dataset 6, frames 2 and 3 from dataset 8, and
% frame 2 from datasets 9 and 10.  Note that data set 6 was not used in
% alaysis but was kept as a reference for 'bad data.' Additionally, the data
% naming must be in ascending numerical order. Discarding all frames avoids
% having to comment out the line and renumber remaining datasets. 
SiMPullUtils.formatFilenames(start_DataDir, ImSZ);

%% ---------- Fit Single Molecules in Each Channel ----------
% Fit the localizations in each channel (using fitblob), then overlay all
% possible channel pairs (and 3 at once if all 3 channels were specified) to
% find all possible colocalizations (using findOverlaps).  The results,
% results.mat and results.xls, are saved for later use in DataDir, which is
% where the images (or tracked files) close allare located.

% Define parameters for each channel.  These will be assigned to the smite SMF
% (Single Molecule Fitting) structure when the images are processed to find
% localization coordinates.  Most of these parameters are set to default values
% as indicated by the comments above their definitions, however, few of these
% have needed to be adapted.

% Pixallow is the number of pixels allowed to be shifted between overlapping
% blobs to consider them as being colocalized.  Note that the blobs are square.
pixallow = 1;

% Various SMITE Single Molecule Fitting (SMF) parameters:

% Overlap of boxes allowed (Pixels)(Default=2)
   smf.BoxFinding_BoxOverlap = [2, 2, 2];
% Linear box size for fitting (Pixels)(Default=7)
   smf.BoxFinding_BoxSize = [5, 5, 5];
% Minimum number of photons from emitter (Default=200)
   smf.BoxFinding_MinPhotons = [75, 75, 75]; 
% Camera Gain, scalar or image (Default=1) ADU / e-
   smf.Data_CameraGain = [17.9, 17.9, 17.9];
% Camera Offset, scalar or image (Default=0) ADU
   smf.Data_CameraOffset = [112.6, 112.6, 112.6];
% See fit class for options: XYNB, XYNBS, XYNBSXSY, XYZNB (Default='XYNB')
  smf.Fitting_FitType = {'XYNBS', 'XYNBS', 'XYNBS'};
% Initial or fixed Sigma of 2D Gaussian PSF Model (Pixels)(Default=1)
   smf.Fitting_PSFSigma = [1, 1, 1];
   % Comment out the above line and uncomment the one below if not using a
   % sigma estimating fit type such as XYNB.
   %smf.Fitting_PSFSigma = [1.2, 1.2, 1.2];
% Maximum allowed precision in x,y (Pixels)(Default=.2)
   smf.Thresholding_MaxXY_SE  = [0.2, 0.2, 0.2];
% Minimum accepted photons from fit (Default=100)
   smf.Thresholding_MinPhotons = [475, 75, 25];
% Minimum accepted p-value from fit (Default=.01)
   smf.Thresholding_MinPValue = [0.01, 0.01, 0.000001];

% Print/save extra information if > 0, such as fit box plots.
Verbose = 1;
SiMPullUtils.overlayFitting(start_DataDir, ImSZ, Indx, Channels,  ...
                            ChannelLabel, pixallow, smf, Verbose, ...
                            UseRegistration);

%% ---------- OPTIONAL: Optimize Minimum Photon Threshold ----------
% To get optimal thresholding after the fitting step, meaning the low signal 
% background autofluorescence are excluded, optimization using this section may
% be needed.  This is particularly important with the 488 data in channel 1, as
% there will be remaining autofluorescence that will skew the final number of
% fits.  We optimize for this by running the fitting and thresholding on the
% blank/background only data using a low minimum photon threshold.  Starting
% with Thresholding_MinPhotons = 0 or 10 worked well for testing.  Next,
% run the section a second time, this time plotting histograms for the
% experimental data. The demo data has been optimized in channel 1, where 
% Thresholding_MinPhoton = 475 exludes most of the background fits.
% 
% NOTE: Setting smf.BoxFinding_MinPhotons = [0, 0, 0] or some other low value
% provides a larger sameple size to see the histogram distribution.  
% 
% This process may need to be iterated a few times.  When a final value is
% decided upon, for which most of the erroneous fits have been eliminated,
% change the corresponding value in the Fit single molecules section as well.
% The output will save a png of the histograms in the Figures folder.  You
% should rename the file produced to distinguish between backgrouund and data
% if you want to save both histograms.

k = 1;   % channel number
smf.Thresholding_MinPhotons = [0, 0, 0];
SiMPullUtils.findMinPhotonThreshold(start_DataDir, k, ImSZ, Indx, Channels, ...
                                    ChannelLabel, smf, Verbose);

%% ---------- OPTIONAL: Optimize Maximum X/Y Standard Error Threshold ----------
% To get optimal thresholding after the fitting step, meaning the low signal 
% background autofluorescence are excluded, optimization using this section may
% be needed.  This is particularly important with the 488 data in channel 1, as
% there will be remaining autofluorescence that will skew the final number of
% fits.  We optimize for this by running the fitting and thresholding on the
% blank/background only data using a high maximum XY_SE threshold.  Starting
% with Thresholding_MaxXY_SE = 0.5 or 1 worked well for testing.  Next,
% run the section a second time, this time plotting histograms for the
% experimental data. The demo data has been optimized in channel 1, where 
% Thresholding_MaxXY_SE = 0.2 exludes most of the background fits.
% 
% NOTE: Setting smf.BoxFinding_MaxXY_SE = [0.5, 0.5, 0.5] or some other high
% value provides a larger sameple size to see the histogram distribution.  
% 
% This process may need to be iterated a few times.  When a final value is
% decided upon, for which most of the erroneous fits have been eliminated,
% change the corresponding value in the Fit single molecules section as well.
% The output will save a png of the histograms in the Figures folder.  You
% should rename the file produced to distinguish between backgrouund and data
% if you want to save both histograms.

k = 1;   % channel number
smf.Thresholding_MaxXY_SE = [0.5, 0.5, 0.5];
SiMPullUtils.findMaxXY_SEThreshold(start_DataDir, k, ImSZ, Indx, Channels, ...
                                   ChannelLabel, smf, Verbose);

%% ---------- Calculate Percentage of GFP fits Positive for FR Signal ----------
% This section determines the 488 background correction and applies it to the
% data in results.mat to calculate the percentage of phosphorylated receptors.
% Before starting, rename results.mat to results_LABEL.mat, where LABEL is a
% date or condition or other identifying characteristic.  File names for
% corresponding conditions are defined in the file 'conditionsCorrectedBG.m'
% that is found in the same directory as results_LABEL.mat.  This file needs to
% be updated based on experimental condition names.  The final outputs are
% saved in the subdirectory Corrected and are as follows:
%    graph_LABEL_TAG.fig          % overlap bar graphs with SEM error bars
%    graph_LABEL_TAG.png          as above, but saved as a .png
%    resultsStats_LABEL_TAG.xls   Excel spreadsheet of all fits and overlaps
%    table_LABEL_TAG.mat          temporary file used to create the Excel file
% where TAG is a condition defined in 'conditionsCorrectedBG.m'.
%
% Note, here background corrections are applied: NGFP = (NLOC - NBG)*SR where
% NGFP is the corrected number of surface GFP receptors,
% NLOC is the total number of emitters localized,
% NBG  is the expected number of background emitters in the area imaged, and
% SR (surface ratio) is the fraction of receptors located at the cell
% surface. This is also defined in 'conditionsCorrectedBG.m'
% In the demo data, pervanadate (PV) treatment is used to yield maximimum 
% phosphorylation throughout the cell, therefore, surface recpeptor correction 
% is not necessary in PV treated controls and the SR was set to 1.
SiMPullUtils.statsCorrectedBG(start_DataDir, Channels);

%% ---------- OPTIONAL: Summary Plots ----------
% Produce various summary plots for the SiMPull analysis.  These include bar
% graphs of the number of fits in each channel and the percentage overlap of
% the overlayed fits.  The plots will be put in DataDir/Figures.  These plots
% can be used to provide a quick overview of the uncorrected comparative
% results, which can come in handy for understanding the final results.  They
% can also serve as a debugging aid if needed.

SiMPullUtils.summaryPlots(start_DataDir, Channels);

%REMOVE?
