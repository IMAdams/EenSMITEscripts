% Specify the path to the folder containing all results
resultsFolder = 'O:\Cell Path\Lidke Lab\IMAdams\iX71 Single Particle Tracking\202408416_alfa_EGFR_EGF_twocolor\Results';
saveFolder = 'O:\Cell Path\Lidke Lab\IMAdams\iX71 Single Particle Tracking\202408416_alfa_EGFR_EGF_twocolor\Results\Plots';

% Get a list of all files in the results folder
fileList = dir(fullfile(resultsFolder, '*_Results.mat'));

% Sort the files to ensure we process them in pairs
fileNames = {fileList.name};
[sortedNames, sortIndex] = sort(fileNames);
fileList = fileList(sortIndex);

% Initialize a cell array to store the paired filenames
pairedFiles = cell(floor(length(fileList)/2), 2);

% Process files in pairs
pairIndex = 1;
for i = 1:2:length(fileList)-1
    % Get the file names for this pair
    file1 = fileList(i).name;
    file2 = fileList(i+1).name;

    % Extract the base name for saving results
    [~, baseName1, ~] = fileparts(file1);
    [~, baseName2, ~] = fileparts(file2);

    % Remove '_Channel1_Results' or '_Channel2_Results' from the base names
    baseName1 = strrep(baseName1, '_Channel1_Results', '');
    baseName2 = strrep(baseName2, '_Channel2_Results', '');

    % Check if the base names match (excluding the channel number)
    if strcmp(baseName1, baseName2)
        pairedFiles{pairIndex, 1} = file1;
        pairedFiles{pairIndex, 2} = file2;
        pairIndex = pairIndex + 1;
    else
        warning('Mismatch in file pairing: %s and %s', file1, file2);
    end
end

% Remove any empty rows from pairedFiles
pairedFiles = pairedFiles(1:pairIndex-1, :);
% Display the paired files
disp('Paired files:');
for i = 1:size(pairedFiles, 1)
    fprintf('Pair %d:\n', i);
    fprintf('  Channel 1: %s\n', pairedFiles{i, 1});
    fprintf('  Channel 2: %s\n', pairedFiles{i, 2});
    fprintf('\n');
end

% Ask user if they want to proceed
userInput = input('Do you want to proceed with plotting? (y/n): ', 's');

if strcmpi(userInput, 'y')
    % Proceed with plotting
    for i = 1:size(pairedFiles, 1)
        try
        file1 = fullfile(resultsFolder, pairedFiles{i, 1});
        file2 = fullfile(resultsFolder, pairedFiles{i, 2});

        % Extract the base name for saving results
        [~, baseName, ~] = fileparts(pairedFiles{i, 1});
        baseName = strrep(baseName, '_Channel1_Results', '');
        
        % create a name for the plot
        displayName = strrep(baseName, '_', '-');

        % Load TR and SMF from both results files
        [TR1, SMF1, TR2, SMF2] = IMA_loadResultsFile(file1, file2);
        
        % Skip if either channel does not have trajectories
        if isempty(TR1) || isempty(TR2)
            warning('Skipping %s: One or both TR structures are empty', displayName);
            continue;  % Skip to the next iteration of the loop
        end

        % Set up MovieParams
        MovieParams3D = smi_vis.GenerateMovies.prepDefaults();

        % Create plot
        PlotFigure3D = figure('Name', ['3D Plot - ' displayName]);
        PlotAxes3D = axes(PlotFigure3D);
        hold(PlotAxes3D, 'on');
          % Adjust plot properties if needed
        title(PlotAxes3D, displayName);

        % Set up MovieMaker
        MovieMaker3D = smi_vis.GenerateMovies(MovieParams3D);
        MovieMaker3D.SMF = SMF1;
        MovieMaker3D.TR = TR1;
        MovieMaker3D.setVitalParams()
        MovieMaker3D.Params.LineOfSite = [-45, 15];
        MovieMaker3D.prepAxes(PlotAxes3D);
        EmptySMD = smi_core.SingleMoleculeData.createSMD();

        % Plot Channel 1 trajectories
        MovieMaker3D.Params.TrajColor = repmat([1, 0, 0], [numel(TR1), 1]); % Red
        h1 = MovieMaker3D.IMA_makeFrame(PlotAxes3D, MovieMaker3D.TR, [], MovieMaker3D.Params, ...
            MovieMaker3D.SMF, EmptySMD, MovieMaker3D.TR(1).NFrames);

        % Add trajectory numbers for Channel 1
        for i = 1:numel(h1)
            if ~isempty(h1(i))
                x = h1(i).XData(1);
                y = h1(i).YData(1);
                z = h1(i).ZData(1);
                text(x, y, z, num2str(i), 'Color', 'r', 'FontSize', 8);
            end
        end


        % Plot Channel 2 trajectories
        MovieMaker3D.Params.TrajColor = repmat([0, 0, 1], [numel(TR2), 1]); % Blue
        MovieMaker3D.SMF = SMF2;
        MovieMaker3D.TR = TR2;
        h2 = MovieMaker3D.IMA_makeFrame(PlotAxes3D, MovieMaker3D.TR, [], MovieMaker3D.Params, ...
            MovieMaker3D.SMF, EmptySMD, MovieMaker3D.TR(1).NFrames);

        % Add trajectory numbers for Channel 2
        for i = 1:numel(h2)
            if ~isempty(h2(i))
                x = h2(i).XData(1);
                y = h2(i).YData(1);
                z = h2(i).ZData(1);
                text(x, y, z, num2str(i), 'Color', 'b', 'FontSize', 8);
            end
        end


        % Highlight overlapping regions
        overlap_threshold = 1; % Distance threshold for overlap in pixels (adjust as needed)
        min_consecutive_frames = 5; % Minimum number of consecutive overlapping frames
        overlap_points = [];

        for l = 1:numel(TR1)
            for j = 1:numel(TR2)
                common_frames = intersect(TR1(l).FrameNum, TR2(j).FrameNum);
                overlap_streak = [];

                for k = 1:length(common_frames)
                    frame = common_frames(k);
                    idx1 = find(TR1(l).FrameNum == frame, 1);
                    idx2 = find(TR2(j).FrameNum == frame, 1);

                    if ~isempty(idx1) && ~isempty(idx2) && ...
                            isfield(TR1(l), 'X') && isfield(TR1(l), 'Y') && ...
                            isfield(TR2(j), 'X') && isfield(TR2(j), 'Y') && ...
                            idx1 <= length(TR1(l).X) && idx1 <= length(TR1(l).Y) && ...
                            idx2 <= length(TR2(j).X) && idx2 <= length(TR2(j).Y)

                        dist = sqrt((TR1(l).X(idx1) - TR2(j).X(idx2))^2 + ...
                            (TR1(l).Y(idx1) - TR2(j).Y(idx2))^2);

                        if dist <= overlap_threshold
                            overlap_streak = [overlap_streak; TR1(l).X(idx1), TR1(l).Y(idx1), frame];
                        else
                            % Check if the previous streak is long enough
                            if size(overlap_streak, 1) >= min_consecutive_frames
                                overlap_points = [overlap_points; overlap_streak];
                            end
                            overlap_streak = [];
                        end
                    end
                end

                % Check the last streak
                if size(overlap_streak, 1) >= min_consecutive_frames
                    overlap_points = [overlap_points; overlap_streak];
                end
            end
        end

        % Plot overlapping points
        if ~isempty(overlap_points)
            h_overlap = plot3(PlotAxes3D, overlap_points(:,1), overlap_points(:,2), overlap_points(:,3), ...
                'go', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
        else
            h_overlap = [];
        end

        % Find a non-empty handle for each channel
        if ~isempty(h1) && ~isempty(h2)
            if ~isempty(h_overlap)
                legend([h1(1), h2(1), h_overlap], 'Channel 1', 'Channel 2', 'Overlap');
            else
                legend([h1(1), h2(1)], 'Channel 1', 'Channel 2');
            end
        else
            warning('Could not create legend due to empty handles');
        end


        % Enable data cursor mode
        dcm = datacursormode(PlotFigure3D);
        set(dcm, 'UpdateFcn', @myupdatefcn);

        % Save the figures
        saveas(PlotFigure3D, fullfile(saveFolder, [baseName '_3D_plot_both_channels.fig']));
        saveas(PlotFigure3D, fullfile(saveFolder, [baseName '_3D_plot_both_channels.png']));

        % Close the figure to free up memory
        close(PlotFigure3D);
    catch ME
        % If an error occurs, display the error message and continue with the next iteration
        warning('Error processing files %s and %s: %s', file1, file2, ME.message);
        end
 
    end
end

if strcmpi(userInput, 'n')
    disp('Plotting cancelled by user.');
end

function txt = myupdatefcn(~, event_obj)
    pos = get(event_obj, 'Position');
    h = get(event_obj, 'Target');
    info = get(h);
    if strcmp(info.Color, [1 0 0]) % Red, Channel 1
        channel = 'Channel 1';
    else % Blue, Channel 2
        channel = 'Channel 2';
    end
    txt = {['X: ', num2str(pos(1))], ...
           ['Y: ', num2str(pos(2))], ...
           ['Frame: ', num2str(pos(3))], ...
           ['Channel: ', channel]};
end

% % Process files in pairs
% for i = 1:2:length(fileList)-1
%     try
%         % Get the file names for this pair
%         file1 = fullfile(resultsFolder, fileList(i).name);
%         file2 = fullfile(resultsFolder, fileList(i+1).name);
%
%         % Extract the base name for saving results
%         [~, baseName, ~] = fileparts(fileList(i).name);
%         shortBaseName = baseName(1:20);
%
%         % Load TR and SMF from both results files
%         [TR1, SMF1, TR2, SMF2] = IMA_loadResultsFile(file1, file2);
%
%         % Set up MovieParams
%         MovieParams3D = smi_vis.GenerateMovies.prepDefaults();
%
%         % Create plot
%         PlotFigure3D = figure('Name', ['3D Plot - ' shortBaseName]);
%         PlotAxes3D = axes(PlotFigure3D);
%         hold(PlotAxes3D, 'on');
%
%         % Set up MovieMaker
%         MovieMaker3D = smi_vis.GenerateMovies(MovieParams3D);
%         MovieMaker3D.SMF = SMF1;
%         MovieMaker3D.TR = TR1;
%         MovieMaker3D.setVitalParams()
%         MovieMaker3D.Params.LineOfSite = [-45, 15];
%         MovieMaker3D.prepAxes(PlotAxes3D);
%         EmptySMD = smi_core.SingleMoleculeData.createSMD();
%
%         % Plot Channel 1 trajectories
%         MovieMaker3D.Params.TrajColor = repmat([1, 0, 0], [numel(TR1), 1]); % Red
%         h1 = MovieMaker3D.IMA_makeFrame(PlotAxes3D, MovieMaker3D.TR, [], MovieMaker3D.Params, ...
%             MovieMaker3D.SMF, EmptySMD, MovieMaker3D.TR(1).NFrames);
%
%         % Plot Channel 2 trajectories
%         MovieMaker3D.Params.TrajColor = repmat([0, 0, 1], [numel(TR2), 1]); % Blue
%         MovieMaker3D.SMF = SMF2;
%         MovieMaker3D.TR = TR2;
%         h2 = MovieMaker3D.IMA_makeFrame(PlotAxes3D, MovieMaker3D.TR, [], MovieMaker3D.Params, ...
%             MovieMaker3D.SMF, EmptySMD, MovieMaker3D.TR(1).NFrames);
%
%         % Highlight overlapping regions
%         overlap_threshold = 1; % Distance threshold for overlap in pixels (adjust as needed)
%         min_consecutive_frames = 5; % Minimum number of consecutive overlapping frames
%         overlap_points = [];
%
%         for l = 1:numel(TR1)
%             for j = 1:numel(TR2)
%                 common_frames = intersect(TR1(l).FrameNum, TR2(j).FrameNum);
%                 overlap_streak = [];
%
%                 for k = 1:length(common_frames)
%                     frame = common_frames(k);
%                     idx1 = find(TR1(i).FrameNum == frame, 1);
%                     idx2 = find(TR2(j).FrameNum == frame, 1);
%
%                     if ~isempty(idx1) && ~isempty(idx2) && ...
%                             isfield(TR1(l), 'X') && isfield(TR1(l), 'Y') && ...
%                             isfield(TR2(j), 'X') && isfield(TR2(j), 'Y') && ...
%                             idx1 <= length(TR1(l).X) && idx1 <= length(TR1(l).Y) && ...
%                             idx2 <= length(TR2(j).X) && idx2 <= length(TR2(j).Y)
%
%                         dist = sqrt((TR1(l).X(idx1) - TR2(j).X(idx2))^2 + ...
%                             (TR1(l).Y(idx1) - TR2(j).Y(idx2))^2);
%
%                         if dist <= overlap_threshold
%                             overlap_streak = [overlap_streak; TR1(l).X(idx1), TR1(l).Y(idx1), frame];
%                         else
%                             % Check if the previous streak is long enough
%                             if size(overlap_streak, 1) >= min_consecutive_frames
%                                 overlap_points = [overlap_points; overlap_streak];
%                             end
%                             overlap_streak = [];
%                         end
%                     end
%                 end
%
%                 % Check the last streak
%                 if size(overlap_streak, 1) >= min_consecutive_frames
%                     overlap_points = [overlap_points; overlap_streak];
%                 end
%             end
%         end
%
%         % Plot overlapping points
%         if ~isempty(overlap_points)
%             h_overlap = plot3(PlotAxes3D, overlap_points(:,1), overlap_points(:,2), overlap_points(:,3), ...
%                 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
%         else
%             h_overlap = [];
%         end
%
%
%         % Add legends
%         % Find a non-empty handle for each channel
%         if ~isempty(h1) && ~isempty(h2)
%             if ~isempty(h_overlap)
%                 legend([h1(1), h2(1), h_overlap], 'Channel 1', 'Channel 2', 'Overlap');
%             else
%                 legend([h1(1), h2(1)], 'Channel 1', 'Channel 2');
%             end
%         else
%             warning('Could not create legend due to empty handles');
%         end
%         % % Add legends
%         % legend(PlotAxes3D, 'Channel 1', 'Channel 2');
%
%         % Adjust plot properties if needed
%         title(PlotAxes3D, baseName);
%
%         % Save the figures
%         saveas(PlotFigure3D, fullfile(saveFolder, [baseName '_3D_plot_both_channels.fig']));
%         saveas(PlotFigure3D, fullfile(saveFolder, [baseName '_3D_plot_both_channels.png']));
%
%         % Close the figure to free up memory
%         close(PlotFigure3D);
%     catch ME
%         % If an error occurs, display the error message and continue with the next iteration
%         warning('Error processing files %s and %s: %s', fileList(i).name, fileList(i+1).name, ME.message);
%     end
% end
