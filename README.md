# densityScatterChart

[![View densityScatterChart on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/95828-densityscatterchart)

Version: 1.0

Density Scatter Chart is a scatter where the color (and/or transparency) of the markers indicates the density of points.

![Example densityScatterChart](/exampleDSC.png)

How to use:
```
x=randn(1000,1);
y=randn(1000,1)
densityScatterChart(x,y);   % create a chart where color varies with density

%% You can set properties using name value pairs
% Create a chart where transparency varies with density:
densityScatterChart(x, y, 'UseColor', false, 'UseAlpha', true);

% Specify a title:
densityScatterChart(x, y, "Title", "My density scatter chart");

%% You can also set properties after making the chart:
d=densityScatterChart(x, y);

% Make a steeper density view
d.DensityExponent = 2;

% Use alpha, but make it a subtle effect by using a small range:
d.AlphaRange = [.2 .8];
```

