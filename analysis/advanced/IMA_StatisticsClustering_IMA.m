% Perform cluster analysis for comparison of experimental conditions.
% NOTE: sma-cluster and PlotSpread should be on your MATLAB path.
%
% Work flow: (note results are put in the subdirectory 'Analysis')
%   Set important parameters:
%            ROI_sizes, Pixel2nm, start_datadir, oneROI
%      ALWAYS run this first MATLAB section before running the later sections
%      so that necessary parameters are set.
%   Define the ROIs:
%            _ResultsStruct.mat (or _Results.mat) -> _ROIs.mat
%      (each image [ResultsStruct or Results file] produces a separate ROIs
%       file containing all the ROIs for that image; the idea is to do all the
%       ROI selection early on so it need not be repeated---each image can have
%       multiple ROIs; note that this deals with files in a single directory.
%       Many time lab members select the ROIs from the SR _Results.mat files,
%       To transfer these to BaGoL ROIs (necessary, as the ROI files are
%       self-contained, meaning they contain the actual coordinates of the
%       localizations in each ROI), it is necessary to transfer the ROI regions
%       to BaGoL localizations.  See below.)
%   POSSIBLY, define the BaGoL ROIs from the previous ROIs and BaGoL results
%            _ROIs.mat + _BaGoL_ResultsStruct.mat -> _BaGoL_ROIs.mat
%      (the two sets of files should be in corresponding order; the BaGoL
%       coordinates replace the coordinates in the original ROI files)
%   POSSIBLY, define the BaGoL ROIs from the previous ROIs and BaGoL results
%            _ROIs.mat + MAPN_*.mat -> _BaGoL_ROIs.mat
%      (the two sets of files should be in corresponding order; the BaGoL
%       coordinates replace the coordinates in the original ROI files)
%   Statistics for a single condition:
%            _ROIs.mat -> ALL_*
%      (compute statistics for a single condition consisting of a number of
%       images all located in the same directory.
%      Need to define algorithm_range, minPts_range (N) and E_range (epsilon).)
%   Combine results from one or more conditions for further processing:
%            multiple ALL_*_results.mat -> [condition]_results.mat
%      (multiple instantiations of a single condition located in multiple
%       directories are collected into a combined instance in a specified
%       directory.
%      Typically easier to do this renaming by hand.  See Combined statistics
%      below for more details on filename expectations.  A typical example of a
%      renamed ALL_*_results file is:
%         RGY#DBSCAN_N=3,E=50
%      )
%   ALTERNATIVELY, copy the files manually and rename them here:
%            multiple ALL_*_results.mat -> [condition]_results.mat
%   Combined statistics for two or more conditions:
%            multiple [condition]_results_mat -> 
%      (multiple conditions are plotted on the same graph and placed in a
%       subdirectory for this particular analysis.  The filenames are assumed
%       to have the structure:
%          experimentalConditions#analysisConditions_results.mat
%       where experimentalConditions contains no #)
%   ALTERNATIVELY, combined statistics for multiple conditions and experiments
%      (as above, where line colors and types for CDF2 plots are set by hand)
%
% Avoid spaces in filenames and conditions!!!
%% ---------- Set important parameters

doHopkins = true;   % Hopkins' test can be time consuming for dense ROIs
% Set to false for very dense ROIs to avoid crashes due to lack of memory.
doSigmaActual = false;

ROI_sizes = [2000, 2000];   % [delta_x, delta_y] (nm)
A_ROI = prod(ROI_sizes);    % ROI area (nm^2)
%Pixel2nm = 16000/150;       % conversion factor from pixels to nm
%Pixel2nm = 108.018;         % pixels to nm [TIRF]
Pixel2nm = 65.0;            % pixels to nm [sequential]
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
start_datadir =  'C:\Users\imadams\Documents\smite workspace\20250822_sorted_CHOHAEGFR\Export_data\Export_data_test_fits';

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

%% ---------- Possibly, define BaGoL ROIs from previous ROIs and BaGoL results
%            (MF BaGoL Results or BaGoL MAPN files)
% NOTE: the two sets of files should be in corresponding order.
% Previously defined ROIs for a series of images.

MAPNfile = false;

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

%% ---------- Possibly, combine individually processed BaGoL ROIs into a single
%           _ROIs.mat file using the _ROIs.mat file that was used to define
%            the ROIs originally from the SR data
[pathnameR, filesR] = smi_helpers.selectFiles(start_datadir, ...
                         '_ROIs.mat file', '*_ROIs.mat');
% BaGoL MAPN or Results/ResultStruct files.
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

%% ---------- Possibly, filter out some ROIs
% Interactive input.
[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
                       '_ROIs.mat files', '*_ROIs.mat');
[n_ROIs, RoI] = CI.filterROIs(pathname, files, filter);
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
E_range = [15];   % (nm)
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

%% ---------- Combine results from 1 or more conditions for further processing

analysis_dir = uigetdir(start_datadir, 'Directory for combined analyses');
if analysis_dir == 0
   error('No directory specified!');
end
Files = uipickfiles('FilterSpec', start_datadir, 'REFilter', ...
                    '.*_results.mat', 'Prompt', '*_results.mat files');
if ~iscell(Files) & Files == 0
   error('No files specified!');
end
[~, in_file, ~] = fileparts(Files{1});
out_file = regexprep(in_file, '^ALL_', '');
answer = inputdlg('Combined filename', '', [1, 50], {out_file});
if numel(answer) == 0
   error('No Filename specified!');
end
out_file = [answer{1}, '.mat'];

CI.combineResults(Files, analysis_dir, out_file);


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
SC.PlotDo = 'nCSx';
% Red mean, green median (2 only mean, 3 only median) for PlotSpread plots.
SC.ShowMM = 1;
% Options for CDF2 plots are: 'plot', 'semilogx', 'semilogy', 'loglog'.
SC.LinLog = 'semilogx';

[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
                       '_results.mat files', '*_results.mat');
answer = inputdlg('Output directory identifier:');
base_name = answer{1};

CI.combinedStatistics1(SC, pathname, files, base_name, A_ROI, doHopkins);

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
files = fliplr(files); % re
answer = inputdlg('Output directory identifier:');
base_name = answer{1};


CI.combinedStatistics1(SC, pathname, files, base_name, ...
                       A_ROI, doHopkins);

