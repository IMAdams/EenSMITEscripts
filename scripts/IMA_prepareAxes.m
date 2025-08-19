function IMA_prepareAxes(PlotAxes, Params, Is3D)
    hold(PlotAxes, 'on')
    PlotAxes.ActivePositionProperty = 'position';
    PlotAxes.DataAspectRatio = [1, 1, 1];
    PlotAxes.YDir = 'reverse';
    PlotAxes.XLimMode = 'manual';
    PlotAxes.YLimMode = 'manual';
    
    if ~isempty(Params.XPixels)
        PlotAxes.XLim = Params.XPixels + [0, 1];
    end
    if ~isempty(Params.YPixels)
        PlotAxes.YLim = Params.YPixels + [0, 1];
    end
    
    if Is3D
        PlotAxes.ZLimMode = 'manual';
        if ~isempty(Params.ZFrames)
            PlotAxes.ZLim = Params.ZFrames;
        end
        view(PlotAxes, Params.LineOfSite(1, :));
    end
    
    colormap(PlotAxes, 'gray')
    
    % Add labels
    xlabel(PlotAxes, 'X (\mu m)', 'Interpreter', 'latex')
    ylabel(PlotAxes, 'Y (\mu m)', 'Interpreter', 'latex')
    if Is3D
        zlabel(PlotAxes, 'Time (frames)', 'Interpreter', 'latex')
    end
end