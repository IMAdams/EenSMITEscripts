% A simple scipt to invoke the plotROIDriver (which invokes plotROI) in various
% ways.

close all

%plotROIDriver(PixelSize, options, start_datadir, SaveDir)

% One BaGoL MAPN image exhibiting multiple ROIs.
%options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};
%options = {'MAPN', 'GaussSEConst', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};

% One BaGoL MAPN image per ROI.
%options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'ROIImages', 'NoSave'};

%options = {'MAPN', 'Gaussian', 'Cluster', 'ROIImages'};
%options = {'MAPN', 'GaussSEConst', 'Cluster', 'ROIImages'};
%options = {'MAPN', 'Dot', 'Cluster', 'ROIImages'};
%options = {'MAPN', 'Circle', 'Cluster', 'ROIImages'};

% One BaGoL MAPN ROI per MAPN Results file.
options = {'MAPN', 'MAPNResultsROI', 'GaussSEConst', 'Cluster'};

% One SMITE SR image exhibiting multiple ROIs.
%options = {'SR', 'Gaussian', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};
%options = {'SR', 'Dot', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};

% One SMITE SR image per ROI.
%options = {'SR', 'Gaussian', 'Boundary', 'Cluster', 'ROIImages', 'NoSave'};
%options = {'SR', 'Gaussian', 'Cluster', 'ROIImages'};

%start_DataDir = '/mnt/nas/cellpath/Genmab/Data/';
start_DataDir = 'O:\Cell Path\Lidke Lab\Genmab\EB Notebook\Data\dSTORM\20250603_A431_SaturatingIgG10min+C1q\Analyses\COMBINED\E430G_C1q_250603\BaGoL_ROIs\DBSCAN_N=3,E=15';

%SaveDir = fullfile(start_DataDir, 'Disc', 'Resting');
%SaveDir = fullfile(start_DataDir, 'Disc', 'DF3');
SaveDir = fullfile(start_DataDir, 'Disc', 'A431_E430G_C1q');

% Cells to produce plots for using an absolute numbering over the cells that
% were selected,
IncludeCell = [];
%IncludeCell = [2]; % Resting Cell 2 ROI 7
%IncludeCell = [3]; % DF3 Cell 1 ROI 11
%IncludeCell = []; % DNP-BSA_10nM Cell 3 ROI 8

%IncludeCell = [13, 25]; % DF3
%IncludeCell = [11, 20]; % DNP-BSA
%IncludeCell = [9, 116]; % DNP-PEG-BSA
%IncludeCell = [1, 8]; % DNP-BSA_10nM
%IncludeCell = [1, 7]; % DNP-BSA_2150pM
%IncludeCell = [2]; % Resting

PixelSize = 97.8;   % nm/pixel

% Plot dot, Gaussian or circle images of SMD/MAPN coordinates per ROI per Cell.
% In other words, this function can plot cluster boundaries overlaid on dot,
% Gaussian or circle plots of SR or BaGoL results.  The driver collects the
% paths and files needed, while plotROI produces the specified plots from these
% and a few other pieces of information for each ROI in each Cell, putting the
% results in SaveDir.  WARNING: Not ALL combinations of options work together!
%
% INPUTS:
%    PixelSize       plxel length (nm)
%    options         cell array of strings that are selected from the following
%                    options [DEFAULT = {'BaGoL', 'Gaussian', 'Cluster'}]:
%       SR              SR Results file
%       BaGoL           BaGoL Results file (BGL.SMD)
%       MAPN            BaGoL MAPN file
%       Dot             Dot plot
%       Gaussian        Gaussian plot
%       GaussSEConst    Gaussian plot with constant X/Y_SE
%       Circle          Circle plot (BaGoL Results file: BGL.SMD + BGL.MAPN)
%       Boundary        Include ROI boundaries
%       Cluster         Include ROI clusters
%       ROIImages       produce one image per ROI
%       OneImage        produce one image displaying all ROIs
%       NoSave          do not save outputs
%    start_datadir   starting data directory from which SaveDir will be made a
%                    subdirectory of by default if SaveDir not provided
%                    [DEFAULT = '.']
%    SaveDir         directory in which to save images
%                    [DEFAULT = fullfile(start_datadir, 'ROIClusterAnalysis')]
%    IncludeCell     Cells to produce plots for [DEFAULT: 1 : numel(filesB)]
%                    using the numbering in the list filesB

CI = smi_cluster.ClusterInterface;
CI.plotROIDriver(PixelSize, options, start_DataDir, SaveDir, IncludeCell);
