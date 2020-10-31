function [ ] = AxesSettings( ax, FontSize, LineWidth, AspectRatio )
%% --------------------------------------------------------------------- %%
if( nargin < 1 )
    ax = findobj( gca );
end
if( nargin < 2 )
    FontSize = 11;
end
if( nargin < 3 )
    LineWidth = 1;
end
%% --------------------------------------------------------------------- %%
ax.LineWidth = LineWidth;
ax.FontSize =  FontSize;
ax.FontName = 'cmr10';
ax.FontWeight = 'normal';
ax.FontAngle = 'normal';
ax.TickLabelInterpreter = 'latex';
ax.TickDir = 'in';
ax.XDir = 'normal';
ax.YDir = 'normal';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.XScale = 'linear';
ax.YScale = 'linear';
ax.Color = [ 1, 1, 1 ];
if( nargin == 4 )
    ax.PlotBoxAspectRatio = [ 1, AspectRatio, 1 ];
else
    ExpandAxesToMinimizeWhiteSpace( ax );
end
% ax.XTick; ax.YTick; ax.XTickLabel; ax.YTickLabe; ax.Units; ax.Position;
%% --------------------------------------------------------------------- %%
end
