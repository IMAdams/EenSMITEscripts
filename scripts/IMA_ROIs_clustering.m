%% ---------- Set important parameters

doHopkins = true;   % Hopkins' test can be time consuming for dense ROIs
% Set to false for very dense ROIs to avoid crashes due to lack of memory.
doSigmaActual = false;

ROI_sizes = [2000, 2000];   % [delta_x, delta_y] (nm)
A_ROI = prod(ROI_sizes);    % ROI area (nm^2)
%Pixel2nm = 16000/150;       % conversion factor from pixels to nm
%Pixel2nm = 108.018;         % pixels to nm [TIRF]
Pixel2nm = 97.8;            % pixels to nm [sequential]
Pixel2nmGlobal = Pixel2nm;
CI = smi_cluster.ClusterInterface();

RT = smi_helpers.ROITools();
% GaussIm = true indicates that gaussianImage will be used for the ROI
%    selection display.
% OriginLLvsUL = true says to use lower left origin coordinates rather than
%    upper left origin coordinates in the ROI selection display.
% If GaussIm is true, make OriginLLvsUL false for consistency.  The default is:
% GaussIm = false and OriginLLvsUL = true.  Currently, GaussIm only works for
% single labeled molecules.
RT.GaussIm = true;    RT.OriginLLvsUL = false;
% GaussIm = true;
% OriginLLvsUL = false;
% RT.GaussIm = false;   RT.OriginLLvsUL = true;
RT.SRzoom = 4;              % zoom factor for gaussianImage
RT.ROI_sizes = ROI_sizes;   
RT.Pixel2nm = Pixel2nm;


% Often, for BaGoL analyses, it is simpler to use a single, large, encompassing
% ROI rather than a series of small ROIs.
oneROI = false;
% oneROI = true;
if oneROI
   ROI_sizes = [256, 256] * Pixel2nm;   % (nm)
   A_ROI = prod(ROI_sizes);
end

% Select the files starting from start_datadir.
start_datadir =  'C:\Users\imadams\Documents\smite workspace\20251024_hela\Bagol_1124';

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
% E_range = [10, 15, 20, 25];   % (nm)
E_range = [18]
% Minimum number of points in a cluster.
%minPts_range = [3, 6, 10];
minPts_range = [3];
% Clustering algorithm.
%algorithm_range = {'DBSCAN', 'Hierarchical', 'Voronoi'};
algorithm_range = {'DBSCAN'};
% Ratio of local density / overall density for Voronoi clustering.
Alpha = 2;

CI.singleCondition(pathname, files, algorithm_range, E_range, minPts_range, ...
                   Pixel2nm, base_name, A_ROI, doHopkins, doSigmaActual, Alpha);

%% ---------- Combined statistics for one or more conditions

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
SC.PlotDo = 'fnpcx';
% Red mean, green median (2 only mean, 3 only median) for PlotSpread plots.
SC.ShowMM = 1;
% Options for CDF2 plots are: 'plot', 'semilogx', 'semilogy', 'loglog'.
SC.LinLog = 'semilogx';

[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
                       '_results.mat files', '*_results.mat');
answer = inputdlg('Output directory identifier:');
base_name = answer{1};

CI.combinedStatistics1(SC, pathname, files, base_name, A_ROI, doHopkins);