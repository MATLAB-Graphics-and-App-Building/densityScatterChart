classdef densityScatterChart < matlab.graphics.chartcontainer.ChartContainer ...
        & matlab.graphics.chartcontainer.mixin.Colorbar
% densityScatterChart - Create a scatter chart that indicates density
% with color or alpha.
%
%     DENSITYSCATTERCHART(x,y) - creates a scatter chart with filled circles
%     where color indicates the density of points. Specify x and y as
%     numeric vectors with matching length.
%
%     DENSITYSCATTERCHART(target,...) - creates the chart in target, for
%     instance a figure or tiled chart layout.
%
%     DENSITYSCATTERCHART(...,Name,Value) - sets the density scatter chart
%     properties using one or more name-value pair arguments. Name/Value
%     pairs are described below.
%
%     dsc = DENSITYSCATTERCHART(...) returns the DensityScatterChart
%     object. Use dsc to set properties on the chart after creating it.
%
%     ---------------------------------------------------------------------
%     DENSITYSCATTERCHART properties description:
%
%     PropertyName default (if relevant) : purpose
%
%     ALim [0 1] : The data limits of alpha used by the chart.
%     ALimMode 'auto' : When 'auto', ALim will match the range of
%                        densities in the chart
%     AlphaRange [0.1 1] : The range of alpha values that ALim will map on
%                           to. 0 inidcates a fully transparent marker, 1
%                           indicates a fully opaque marker.
%     CLim [0 1] : The data limits of color used by the chart.
%     CLimMode 'auto' : When 'auto', CLim will match the range of
%                       densities in the chart.
%     ColorbarVisible on : Whether or not the colorbar is visible.
%     Colormap : The colormap used by the chart. Defaults to the default
%                colormap used by MATLAB
%     DensityExponent 1 : An exponent applied to density values. Use a
%                         value greater than 1 for a 'steeper' density
%                         display, and a value less than 1 for a 'flatter'
%                         density display.
%     DensityMethod "histcounts" : The method used to compute density.
%                                    - "histcounts" will bin the data in
%                                      rectangles (using the histcounts2
%                                      function) and density will be the
%                                      number of points in each bin.
%                                    - "ksdensity" requires the Statistics
%                                      and Machine Learning toolbox, and will
%                                      calculate the kernel density (using
%                                      the ksdensity function with the
%                                      default arguments).
%                                    - You may also specify a custom function
%                                      that takes two arguments (the x and y
%                                      values) and returns a density of matching
%                                      length, e.g.: @(x,y)x^2+y^2;
%     UseAlpha off : Whether or not to vary alpha with density.
%     UseColor on : Whether or not to vary color with density.
%     XData: [1×0 double]
%     YData: [1×0 double]
%     Title "" : A title for the chart.
%
%     Additional properties (see documentation for axes)
%       Position, InnerPosition, OuterPosition, PositionConstraint, Parent,
%       Units, Visible, XLim, XLimMode, YLim, YLimMode

% Copyright 2021 The MathWorks, Inc.

    % Public interface:
    properties
        XData (1,:) double = []
        YData (1,:) double = []
        Colormap (:,3) double {mustBeNonempty, mustBeInRange(Colormap,0,1)} = get(groot, 'factoryFigureColormap')
        UseColor (1,1) matlab.lang.OnOffSwitchState = matlab.lang.OnOffSwitchState.on
        UseAlpha (1,1) matlab.lang.OnOffSwitchState = matlab.lang.OnOffSwitchState.off
        AlphaRange (1,2) double {mustBeLimits} = [.1 1]

        DensityExponent (1,1) double {mustBePositive} = 1;
        DensityMethod {mustBeDensityMethod} = "histcounts"

        Title (:,1) string = ""
    end

    % DataStorage, used for save/load
    properties(Access = protected)
        DataStorage struct
    end

    % Limits properties that 'live' on the axes
    properties(Dependent)
        XLim
        YLim
        CLim
        ALim
        XLimMode
        YLimMode
        CLimMode
        ALimMode
    end

    properties(Access = private, Transient, NonCopyable)
        Scat (1,1) matlab.graphics.chart.primitive.Scatter

        DataNeedsUpdate (1,1) logical = true
        DensityNeedsUpdate (1,1) logical = true
    end

    % Chart constructor
    methods
        function obj = densityScatterChart(varargin)
            args = varargin;

            % Check if the first input argument is a graphics object to use as parent.
            leadingArgs = cell(0);
            if ~isempty(args) && isa(args{1},'matlab.graphics.Graphics')
                % densityScatterChart(parent, ___)
                leadingArgs = args(1);
                args = args(2:end);
            end

            % Check for optional positional arguments.
            if ~isempty(args) && numel(args) >= 2 && ...
                    isnumeric(args{1}) && isnumeric(args{2})
                % densityScatterChart(x, y)
                % densityScatterChart(x, y, Name, Value)
                x = args{1};
                y = args{2};
                leadingArgs = [leadingArgs {'XData', x, 'YData', y}];
                args = args(3:end);

                if ~isvector(x) || ~isvector(y) || numel(x) ~= numel(y)
                    throw(MException('densityScatterChart:XYMismatch', ...
                        'The x and y arguments must be vectors of the same length.'))
                end
            end

            if ~isempty(args) && (mod(numel(args),2) == 1 || (~ischar(args{1}) && ~isstring(args{1})))
                throw(MException('densityScatterChart:InvalidArguments', ...
                    'Invalid arguments. Call densityScatterChart as densityScatterChart(x,y), densityScatterChart(...,Name,Value) or densityScatterChart(parent,...). Specify x and y as numeric values.'))
            end


            % Combine positional arguments with name/value pairs.
            args = [leadingArgs args];

            obj@matlab.graphics.chartcontainer.ChartContainer(args{:});
        end
    end

    % Protected chart setup, update, and disp methods
    methods (Access = protected)
        function setup(obj)
            % Create the scatter object
            obj.Scat = matlab.graphics.chart.primitive.Scatter(...
                'Parent', obj.getAxes, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', 'flat', ...
                'Marker', 'o', ...
                'SizeData', 36, ...
                'AlphaDataMapping','scaled');

            box(obj.getAxes, 'on');

            obj.loadState;
        end
        function update(obj)
            if numel(obj.XData) ~= numel(obj.YData)
                warning('XData must match YData')
                return
            end

            valid = isfinite(obj.XData) & isfinite(obj.YData);
            x = obj.XData(valid);
            y = obj.YData(valid);

            if obj.DataNeedsUpdate
                set(obj.Scat, 'XData', x, 'YData', y);
                obj.DensityNeedsUpdate = true;
                obj.DataNeedsUpdate = false;
            end

            obj.getAxes.Alphamap = linspace(obj.AlphaRange(1), obj.AlphaRange(2),64);
            if obj.DensityNeedsUpdate
                d = obj.getDensityData(x,y) .^ obj.DensityExponent;
                if obj.UseColor
                    % Color the values based on density.
                    obj.Scat.CData = d;
                    obj.ColorbarVisible = 'on';
                else
                    % If not coloring by density, choose a color from the
                    % default axes colororder.
                    if ~isempty(obj.Parent)
                        c=get(obj, 'DefaultAxesColorOrder');
                    else
                        c=get(groot, 'DefaultAxesColorOrder');
                    end
                    obj.Scat.CData = c(1,:);

                    % Workaround for issues with Colorbar mixin
                    try
                        obj.ColorbarVisible = 'off';
                    end
                end
                if obj.UseAlpha
                    % Set the alpha based on density, but flip it!
                    obj.Scat.MarkerFaceAlpha = 'flat';
                    obj.Scat.AlphaData = max(d)-d;
                else
                    obj.Scat.MarkerFaceAlpha = 1;
                end
            end
            title(obj.getAxes, obj.Title);
        end
        function groups = getPropertyGroups(~)
            groups = matlab.mixin.util.PropertyGroup( ...
                {'XData', 'YData', ...
                'UseColor', 'UseAlpha', ...
                'DensityMethod', 'DensityExponent'});
        end

        % The density calculation
        function d = getDensityData(obj,x,y)
            if isempty(obj.XData)
                d = [];
                return
            end
            if numel(x)<3
                d=ones(size(x));
                return
            end

            if strcmpi(obj.DensityMethod,'ksdensity')
                % stats toolbox density via ksdensity
                d = ksdensity([x y], [x y]);

            elseif strcmpi(obj.DensityMethod,'histcounts')
                [n,xedges,yedges] = histcounts2(x, y, 'BinMethod', 'auto');

                xcenters = xedges(1:end-1) + diff(xedges) / 2;
                ycenters = yedges(1:end-1) + diff(yedges) / 2;

                if numel(ycenters) == 1
                    d = interp1(xcenters,n',x);
                elseif numel(xcenters) == 1
                    d = interp1(ycenters,n',y);
                else
                    [xi,yi] = meshgrid(xcenters, ycenters);
                    d = interp2(xi, yi, n', x, y);
                end

            elseif isa(obj.DensityMethod, 'function_handle')
                try
                    d = obj.DensityMethod(x, y);
                catch ME
                    d = ones(size(x));
                    warning("Error in DensityMethod: " + ME.message)
                end
            end
        end
    end



    % Set and Get methods
    methods

        % Data Properties: should mark the DataNeedsUpdate flag when these
        % change
        function set.XData(obj, val)
            obj.XData = val;
            obj.DataNeedsUpdate = true;
        end
        function set.YData(obj, val)
            obj.YData = val;
            obj.DataNeedsUpdate = true;
        end

        % Limits and accompanying mode properties get passed through to
        % axes
        function set.XLim(obj, val)
            obj.getAxes.XLim = val;
        end
        function set.YLim(obj, val)
            obj.getAxes.YLim = val;
        end
        function set.CLim(obj, val)
            obj.getAxes.CLim = val;
        end
        function set.ALim(obj, val)
            obj.getAxes.ALim = val;
        end
        function set.XLimMode(obj, val)
            obj.getAxes.XLimMode = val;
        end
        function set.YLimMode(obj, val)
            obj.getAxes.YLimMode = val;
        end
        function set.CLimMode(obj, val)
            obj.getAxes.CLimMode = val;
        end
        function set.ALimMode(obj, val)
            obj.getAxes.ALimMode = val;
        end

        function val = get.XLim(obj)
            val = obj.getAxes.XLim;
        end
        function val = get.YLim(obj)
            val = obj.getAxes.YLim;
        end
        function val = get.CLim(obj)
            val = obj.getAxes.CLim;
        end
        function val = get.ALim(obj)
            val = obj.getAxes.ALim;
        end
        function val = get.XLimMode(obj)
            val = obj.getAxes.XLimMode;
        end
        function val = get.YLimMode(obj)
            val = obj.getAxes.YLimMode;
        end
        function val = get.CLimMode(obj)
            val = obj.getAxes.CLimMode;
        end
        function val = get.ALimMode(obj)
            val = obj.getAxes.ALimMode;
        end

        % datastorage property supports saving and loading the chart
        function data=get.DataStorage(obj)
            % this method is called when the chart is saved or loaded. It
            % stores properties that 'live' on the axes
            isLoading = ~isempty(obj.DataStorage);
            if isLoading
                data = obj.DataStorage;
            else
                data = struct('XLim', [], 'YLim', [], 'CLim', [], 'ALim', []);
                if strcmp(obj.getAxes.XLimMode, 'manual')
                    data.XLim = obj.getAxes.XLim;
                end
                if strcmp(obj.getAxes.YLimMode, 'manual')
                    data.YLim = obj.getAxes.YLim;
                end
                if strcmp(obj.getAxes.CLimMode, 'manual')
                    data.CLim = obj.getAxes.CLim;
                end
                if strcmp(obj.getAxes.ALimMode, 'manual')
                    data.ALim = obj.getAxes.ALim;
                end
            end
        end
    end

    % A custom method is used to support axes properties during load
    methods (Access = protected)
        function loadState(obj)
            data = obj.DataStorage;
            if isempty(data)
                return
            end

            f = fieldnames(data);
            for i = 1:numel(f)
                fn = f{i};
                if ~isempty(data.(fn))
                    obj.getAxes.(fn) = data.(fn);
                end
            end
            obj.DataStorage = [];
        end
    end

    % Support for MATLAB convenience functions
    methods
        function out = colormap(obj, varargin)
            o = colormap(obj.getAxes, varargin{:});
            if nargout == 1
                out = o;
            end
        end
        function out = xlim(obj, varargin)
            if nargout == 1
                out = xlim(obj.getAxes, varargin{:});
            else
                xlim(obj.getAxes, varargin{:});
            end
        end
        function out = ylim(obj, varargin)
            if nargout == 1
                out = ylim(obj.getAxes, varargin{:});
            else
                ylim(obj.getAxes, varargin{:});
            end
        end
        function out = caxis(obj,varargin)
            if nargout == 1
                out = caxis(obj.getAxes, varargin{:});
            else
                caxis(obj.getAxes, varargin{:});
            end
        end
    end

end

% Property Validators
function mustBeLimits(a)
if numel(a) ~= 2 || a(2) <= a(1) || any(a < 0) || any(a > 1)
    throwAsCaller(MException('densityScatterChart:InvalidLimits', 'Specify alpha range as two increasing values between 0 and 1.'))
end
end

% Property Validators
function mustBeDensityMethod(a)
if isa(a, 'function_handle')
    if nargin(a) == 2
        return
    else
        throwAsCaller(MException('densityScatterChart:InvalidDensityFunc', ...
            'When specifying ''DensityMethod'' as a function handle, the function must accept two arguments.'))
    end
end

if ischar(a) || isstring(a)
    if strcmpi(a, 'histcounts')
        return
    end
    if strcmpi(a, 'ksdensity')
        if license('test', 'statistics_toolbox')
            return
        else
            throwAsCaller(MException('densityScatterChart:InvalidDensityStatsToolbox', ...
                'Statistics toolbox must be installed to used the ksdensity option.'))
        end
    end
end

throwAsCaller(MException('densityScatterChart:InvalidDensity', ...
    '''DensityMethod'' must be a function handle, or the keywords ''histcounts'' or ''kdensity''.'))
end
