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
RT.GaussIm = true; RT.OriginLLvsUL = true;
% RT.GaussIm = true; RT.OriginLLvsUL = false;
RT.SRzoom = 4;         
RT.ROI_sizes = ROI_sizes;   
RT.Pixel2nm = Pixel2nm;  % previously the RT obj had a default of 100.0 nm 
oneROI = false;

start_datadir =  'C:\Users\imadams\Documents\smite workspace\ROI cluster overlays';
%start_datadir =  '/mnt/nas/cellpath/Personal Folders/IMAdams/Overlay Error';

fprintf('Done set parameters.\n');

%% ----------- Define the ROIs
[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
                       '_Results*.mat files', '*_Results.mat');
%files = uipickfiles(FilterSpec, start_datadir, REFilter, ...
%                    '_Results*\.mat', 'Prompt', '_Results*.mat files');
results_dir = fullfile(pathname, 'Analysis');
% Create results_dir if it does not already exist.
if ~isfolder(results_dir)
   mkdir(results_dir);
end
% Define the ROIs for each image.
Pixel2nm = Pixel2nmGlobal;
n_files = numel(files);
for j = 1 : n_files
   short = regexprep(files{j}, '.mat$', '');
   ResultsFile = fullfile(pathname, files{j});
   fprintf('%s ...\n', short);
   if oneROI
      [XY, XY_STD, XYsize] = ...
         RT.import_XY(ResultsFile, Pixel2nm, '');
      n_ROIs = 1;
      RoI{1}.ROI = [0, ROI_sizes(1), 0, ROI_sizes(2)];
      RoI{1}.X = {XY(:, 1)};
      RoI{1}.Y = {XY(:, 2)};
      RoI{1}.X_STD = {XY_STD(:, 1)};
      RoI{1}.Y_STD = {XY_STD(:, 2)};
   else
      [n_ROIs, RoI, XYsize] = ...
         RT.getROI(ResultsFile, files{j});
      saveas(gcf, fullfile(results_dir, sprintf('%s_ROIs.fig', short)));
      saveas(gcf, fullfile(results_dir, sprintf('%s_ROIs.png', short)));
   end
   save(fullfile(results_dir, [short, '_ROIs.mat']), ...
        'ResultsFile', 'Pixel2nm', 'XYsize', 'n_ROIs', 'RoI');
   close
end
fprintf('Done ROIs.\n');

%% ---------- Possibly, combine individually processed BaGoL ROIs into a single
%%            _ROIs.mat file using the _ROIs.mat file that was used to define
%%            the ROIs originally from the SR data

[pathnameR, filesR] = smi_helpers.selectFiles(start_datadir, ...
                         '_ROIs.mat file', '*_ROIs.mat');
MAPNfile=true% BaGoL MAPN or Results/ResultStruct files.
if MAPNfile
   [pathnameB, filesB] = smi_helpers.selectFiles(pathnameR, ...
                            'MAPN*_ROI_*.mat files', 'MAPN_*_ROI_*.mat');
else
   [pathnameB, filesB] = smi_helpers.selectFiles(pathnameR, ...
                            '*_Results_ROI_*.mat files',    ...
                            '*_Results_ROI_*.mat');
end

CI.combineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile, ...
                    keep_numbering);

%% ---------- Statistics for a single condition

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
%E_range = [30, 40, 50];   % (nm)
E_range = [50];   % (nm)
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
%start_DataDir = '/mnt/nas/cellpath/Personal Folders/IMAdams/Overlay Error';
PixelSize = 97.4;
%SaveDir = 'C:\Users\imadams\Documents\smite workspace\ROI cluster overlays\wtRestOverlay';
SaveDir = fullfile(start_DataDir, 'wtRestOverlay');

% Cells to produce plots for using an absolute numbering over the cells that
% were selected,
IncludeCell = [];

CI = smi_cluster.ClusterInterface;
CI.plotROIDriver(PixelSize, options, start_DataDir, SaveDir, IncludeCell);
