classdef IMA_TwoChannelMovieGenerator
    properties
        % Structure of parameters enforced in the generated movie.
        Params struct = struct([])
        % Raw data displayed under trajectories. (YSizexXSizex3xNFrames)
        RawData1
        % Single Molecule Fitting structure, channel 1.
        SMF1 = smi_core.SingleMoleculeFitting;
        % Optional SMD containing points to mark in the movie, channel 1.
        SMD1 = smi_core.SingleMoleculeData.createSMD();
        % Tracking Results structure for the trajectories, channel 1.
        TR1 = smi_core.TrackingResults.createTR();
        % similar for channel 2
        RawData2
        TR2 = smi_core.TrackingResults.createTR();
        SMF = smi_core.SingleMoleculeFitting;
        SMD2 = smi_core.SingleMoleculeData.createSMD();
    end

    properties (SetAccess = 'protected')
        % Graphics handles to trajectory lines displayed in movie frames.
        % NOTE: This is index the same way as obj.TR, i.e., obj.TR(n)
        %       corresponds to obj.LineHandles(n).
        LineHandles

        % Figure containing the GUI.
        GUIFigure

        % Axes in which the movie is prepared if using obj.gui().
        MovieAxes

        % Rescaled/cropped version of property RawData (see rescaleData())
        % NOTE: I've made this a protected property so that we can ensure
        %       some of obj.Params are updated when this is written to
        %       (which I only do in obj.rescaleData())
        ScaledData1 % channel 1
        ScaledData2 % channel 2

        % MATLAB VideoWriter object used to write a movie to a file.
        VideoObject
        % Flag indicating 'ScaledData' is ready to use.
        DataIsPrepped = false;

        % Flag indicating 'MovieAxes' is ready to use (e.g., labeled).
        AxesPrepped = false;
    end

    properties (Hidden)
        DispPlotsOptions = {'DatasetNum', 'FrameNum', 'X', 'Y', 'Z', ...
            'X_SE', 'Y_SE', 'Z_SE', ...
            'Photons', 'Photons_SE', 'Bg', 'Bg_SE', ...
            'PSFSigma', 'PSFSigmaX', 'PSFSigmaY', ...
            'PSFSigma_SE', 'PSFSigmaX_SE', 'PSFSigmaY_SE', ...
            'PValue', 'LogLikelihood', 'ThreshFlag', 'ConnectID'}
        LengthUnitOptions = {'pixels', '$\mu m$'};
        TimeDimensionOptions = {'Frame', 'Time'};
        TimeUnitOptions = {'frames', 's'};
    end

    properties (Hidden, Dependent)
        LengthUnitString
        TimeDimensionString
        TimeUnitString
    end

    properties (Hidden, SetAccess = 'protected')
        % Internally modified version of the user set field 'TR'.
        TRInternal = smi_core.TrackingResults.createTR();
    end


    methods
        function obj = IMA_TwoChannelMovieGenerator(MovieParams)
            % Constructor for IMA_TwoChannelMovieGenerator

            % Set default parameters if not provided
            if nargin < 1 || isempty(MovieParams)
                MovieParams = obj.prepDefaults();
            end

            % Initialize parameters
            obj.Params = MovieParams;

            % Initialize other properties
            obj.SMF1 = smi_core.SingleMoleculeFitting;
            obj.SMD1 = smi_core.SingleMoleculeData.createSMD();
            obj.TR1 = smi_core.TrackingResults.createTR();

            obj.SMF2 = smi_core.SingleMoleculeFitting;
            obj.SMD2 = smi_core.SingleMoleculeData.createSMD();
            obj.TR2 = smi_core.TrackingResults.createTR();

            obj.DataIsPrepped = false;
            obj.AxesPrepped = false;

            % Initialize other properties as needed
            obj.LineHandles = [];
            obj.GUIFigure = [];
            obj.MovieAxes = [];
            obj.ScaledData1 = [];
            obj.ScaledData2 = [];
            obj.VideoObject = [];
        end
        function IMA_generateMovie(obj)
            obj.prepRawData()
            obj.prepAxes()

            obj.IMA_playMovieTwoChannels(obj.MovieAxes, obj.TR, obj.TR2, ...
                obj.ScaledData, obj.ScaledData2, obj.Params, obj.SMF, obj.SMF2, ...
                obj.SMD, obj.SMD2, obj.VideoObject);
        end

        function [LengthUnitString] = get.LengthUnitString(obj)
            LengthUnitString = smi_helpers.arrayMUX(...
                obj.LengthUnitOptions, obj.Params.UnitFlag);
        end

        function [TimeDimensionString] = get.TimeDimensionString(obj)
            TimeDimensionString = smi_helpers.arrayMUX(...
                obj.TimeDimensionOptions, obj.Params.UnitFlag);
        end

        function [TimeUnitString] = get.TimeUnitString(obj)
            TimeUnitString = smi_helpers.arrayMUX(...
                obj.TimeUnitOptions, obj.Params.UnitFlag);
        end

        prepRawData(obj)
        saveMovie(obj, SavePath)
        generateMovie(obj)
        gui(obj)
        IMA_setVitalParams(obj, defaultColorName);
        playMovie(PlotAxes, TR, ScaledData, Params, SMF, SMD, VideoObject)
        saveRawDataMovie(RawData, FilePath, Params, FrameRate);
        IMA_playMovieTwoChannels(PlotAxes, TR1, TR2, ScaledData1, ScaledData2, ...
            Params, SMF1, SMF2, SMD1, SMD2, VideoObject);
        IMA_makeFrameTwoChannels(PlotAxes, TR1, TR2, ScaledData1, ScaledData2, ...
            Params, SMF1, SMF2, SMD1, SMD2, Frame);
        IMA_prepRawdata(obj);
        prepAxes(obj, PlotAxes)

    end

    methods (Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        addAxesTicks(obj, PlotAxes, NTicks, FormatSpec)
    end


    methods (Static)
        prepDefaults();

    end

    methods (Static, Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        [LineHandles] = plotTrajectories(PlotAxes, ...
            Params, TR, FrameRange, Color, varargin);
        [LineHandles] = IMA_makeFrame(PlotAxes, TR, ScaledData, ...   % IMA added method for custom function
            Params, SMF, SMD, Frame);
        [LineHandles] = makeFrame(PlotAxes, TR, ScaledData, ...
            Params, SMF, SMD, Frame);
        [Params] = defineCropROI(TR, Params);
        [RawData] = cropRawData(RawData, Params);
        addTimeStamp(PlotAxes, Frame, FrameRate, Params)
    end

end