function [LineHandles] = IMA_makeFrameTwoChannels(PlotAxes, TR1, TR2, ScaledData1, ScaledData2, ...
    Params, SMF1, SMF2, SMD1, SMD2, Frame)

% Clear the axes
cla(PlotAxes)

% Display the raw data from both channels side by side
totalWidth = size(ScaledData1, 2) + size(ScaledData2, 2);
surface(PlotAxes, [1, totalWidth], Params.YPixels + [0, 1], ...
    repmat(min(Params.ZFrames), [2, 2]), ...
    cat(2, ScaledData1(:,:,:,Frame), ScaledData2(:,:,:,Frame)), ...
    'facecolor', 'texturemap', 'EdgeColor', 'none')

hold(PlotAxes, 'on')

% Plot trajectories for channel 1
LineHandles1 = smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
    Params, TR1, [Frame-Params.MaxTrajLength, Frame], ...
    [1, 0, 0], 'Marker', Params.PlotMarker);  % Red for channel 1

% Plot trajectories for channel 2 (adjust X coordinates)
TR2Adjusted = TR2;
for i = 1:numel(TR2)
    TR2Adjusted(i).X = TR2(i).X + size(ScaledData1, 2);
end
LineHandles2 = smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
    Params, TR2Adjusted, [Frame-Params.MaxTrajLength, Frame], ...
    [0, 0, 1], 'Marker', Params.PlotMarker);  % Blue for channel 2

% Plot SMD points if available (with adjustment for channel 2)
if ~isempty(SMD1)
    smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
        Params, SMD1, [Frame-Params.MaxTrajLength, Frame], ...
        [1, 0.5, 0.5], 'Marker', Params.SMDMarker, 'LineStyle', 'none');
end
if ~isempty(SMD2)
    SMD2Adjusted = SMD2;
    SMD2Adjusted.X = SMD2.X + size(ScaledData1, 2);
    smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
        Params, SMD2Adjusted, [Frame-Params.MaxTrajLength, Frame], ...
        [0.5, 0.5, 1], 'Marker', Params.SMDMarker, 'LineStyle', 'none');
end

% Add timestamp if needed
if Params.AddTimeStamp
    smi_vis.GenerateMovies.addTimeStamp(PlotAxes, Frame, ...
        SMF1.Data.FrameRate, Params);
end

% Combine line handles
LineHandles = [LineHandles1; LineHandles2];
end