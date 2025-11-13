%% Code to replicate overlay registration error. I defined RT.Pixel2nm and that improved things but it is still off. 

% ---------- Set important parameters

doHopkins = true;
doSigmaActual = false;

ROI_sizes = [2000, 2000];   % [delta_x, delta_y] (nm)
A_ROI = prod(ROI_sizes);    % ROI area (nm^2)
Pixel2nm = 97.4;            % pixels to nm [sequential]
Pixel2nmGlobal = Pixel2nm;
CI = smi_cluster.ClusterInterface();
RT = smi_helpers.ROITools();
%RT.GaussIm = false; RT.OriginLLvsUL = true;
RT.GaussIm = true; RT.OriginLLvsUL = false;
RT.SRzoom = 4;         
RT.ROI_sizes = ROI_sizes;   
RT.Pixel2nm = Pixel2nm;  % previously the RT obj had a default of 100.0 nm 
oneROI = false;

start_datadir =  'C:\Users\imadams\Documents\smite workspace\ROI cluster overlays';
%start_datadir =  '/mnt/nas/cellpath/Personal Folders/IMAdams/Overlay Error';

% If true, look for MAPN_*.mat, otherwise *_Results*.mat for BaGoL coordinates.
MAPNfile = true;

% keep_numbering retains the ROI numbering even if there are missing ROIs
% (which will be treated as empty).  See combineBaGoLROIs.
keep_numbering = false;

fprintf('Done set parameters.\n');

%% ----------- Define the ROIs
[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
                       '_Results*.mat files', '*_Results.mat');
CI.defineROIs(pathname, files, Pixel2nmGlobal, RT, oneROI);

%% ---------- Possibly, define BaGoL ROIs from previous ROIs and BaGoL results
%%            (MF BaGoL Results or BaGoL MAPN files)

% NOTE: the two sets of files should be in corresponding order.
% Previously defined ROIs for a series of images.
[pathnameR, filesR] = smi_helpers.selectFiles(start_datadir, ...
                           '_ROIs.mat files', '*_ROIs.mat');
% BaGoL MAPN or Results/ResultStruct files.
if MAPNfile
   [pathnameB, filesB] = smi_helpers.selectFiles(start_datadir, ...
                            'MAPN*.mat files', 'MAPN_*.mat');
else
   [pathnameB, filesB] = smi_helpers.selectFiles(start_datadir, ...
                            '*_Results*.mat files', '*_Results*.mat');
end

CI.defineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile, ...
                   RT.OriginLLvsUL);

%% ---------- Statistics for a single condition

% Run stats in order you want them plotted. make sure line 162 with the
% fliplr() function is working. 

% Main parameters to change: algorithm_range, E_range, minPts_range

%  - - - - - Choose _ROIs.mat files

% Interactive input.
[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
                       '_ROIs.mat files', '*_ROIs.mat');
answer = inputdlg('Name of combined results file:');
condition = answer{1};

% Look for condition/Analysis/*_ROIs.mat or uncomment "files = ..."
% lines below and list files manually.
% Use specified condition.
%condition = 'RGY';
%pathname = fullfile(start_datadir, condition, 'Analysis');
%FILES = dir(fullfile(pathname, '*_Results_ROIs.mat'));
%files = { FILES.name };

% Manual specification.
%files = {
%'Cell_01_Label_01_Results_ROIs.mat'
%};

% base_name is a short identifier describing the analysis.
base_name = condition;

%  - - - - - Set up batch analysis ranges

% Minimum distance between clusters or maximum distance between points in a
% cluster.
E_range = [15, 20, 25, 30, 40, 50];   % (nm)
% E_range = [40];   % (nm)
% Minimum number of points in a cluster.
%minPts_range = [3, 6, 10];
minPts_range = [2];
% Clustering algorithm.
%algorithm_range = {'DBSCAN', 'Hierarchical', 'Voronoi'};
algorithm_range = {'DBSCAN'};
% Ratio of local density / overall density for Voronoi clustering.
Alpha = 2;

CI.singleCondition(pathname, files, algorithm_range, E_range, minPts_range, ...
                   Pixel2nm, base_name, A_ROI, doHopkins, doSigmaActual, Alpha);

%% Generate ROI overlays with clusters

close all

options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};
start_DataDir = 'C:\Users\imadams\Documents\smite workspace\ROI cluster overlays';
% start_DataDir = '/mnt/nas/cellpath/Personal Folders/IMAdams/Overlay Error';
PixelSize = 97.4;
%SaveDir = 'C:\Users\imadams\Documents\smite workspace\ROI cluster overlays\wtRestOverlay';
SaveDir = fullfile(start_DataDir, 'wtRestOverlay');

% Cells to produce plots for using an absolute numbering over the cells that
% were selected,
IncludeCell = [];

CI = smi_cluster.ClusterInterface;
CI.plotROIDriver(PixelSize, options, start_DataDir, SaveDir, IncludeCell);
%% ---------- Combined statistics for multiple conditions and experiments

SC = smi_cluster.StatisticsClustering();
% Make various plots:
%    'f'   frequency
%    'n'   normalized
%    'p'   PDF
%    'c'   CDF
%    'C'   CDF (alternative)
%    's'   PlotSpread
%    'S'   PlotSpread (bars for mean & median)
%    'x'   box
%    'b'   bar
SC.PlotDo = 'CSx';
% Red mean, green median (2 only mean, 3 only median) for PlotSpread plots.
%SC.ShowMM = 1;
% Options for CDF2 plots are: 'plot', 'semilogx', 'semilogy', 'loglog'.
SC.LinLog = 'semilogx';
SC.Ylim = [0.01, 1];
SC.CSV = true;
SC.Font_props(2) = {[50]}; % 
% Use the same color for each experiment under the same condition, different
% colors for different conditions.
%colors = ['b', 'b', 'b', 'b', 'r', 'r', 'r', 'r', 'g', 'g', 'g', 'g']; %GFP
% colors = ['b', '#4a5782', 'm', 'k']; % Cyan: EGFR-WT-REST, Dark blue: EGFR-WT +EGF, Magenta: EGFR-L858R Resting, Black: EGFR-L858R +EGF
%colors = ['b', 'b', 'r', 'r', 'g', 'g']; %GFP_exp1+2
%colors = ['b', 'r', 'g']; %GFP3+4combined
%colors = ['b', 'r']; %EGFR

% Use different line types to distingush same colored lines when desired.
%line_type = {':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.', '-', '--'}; %GFP
%line_type = {'-', '-', '-', '--', '--', '--'}; %CD277+parental
%line_type = {'-', '--', '-', '--', '-', '--'}; %GFP_exp1+2
%line_type = {'-', '-', '-'}; %GFP3+4combined
%line_type = {'-', '-'}; %EGFR


[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
                       '_results.mat files', '*_results.mat');
answer = inputdlg('Output directory identifier:');
files = fliplr(files); % if ROIs set in group order
base_name = answer{1};

CI.combinedStatistics1(SC, pathname, files, base_name, ...
                       A_ROI, doHopkins);