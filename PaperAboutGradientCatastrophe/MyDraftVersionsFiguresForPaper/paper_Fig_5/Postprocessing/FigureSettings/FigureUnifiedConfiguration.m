function [ ] = FigureUnifiedConfiguration( fig, FontSize, LineWidth )
%% --------------------------------------------------------------------- %%
if( nargin < 1 )
    fig = findobj( gcg );
end
if( nargin < 2 )
    FontSize = 11;
end
if( nargin < 3 )
    LineWidth = 1;
end
%% --------------------------------------------------------------------- %%
allaxes  = findall( fig, 'Type', 'axes' );
alltext  = findall( fig, 'Type', 'text' );
alllines = findall( fig, 'Type', 'line' );
set(alllines, 'LineWidth', LineWidth );
set( allaxes, 'FontSize', FontSize, ...
              'FontName', 'cmr10', 'FontWeight', 'normal', 'FontAngle', 'normal' );
set( alltext, 'Interpreter', 'latex', 'FontSize', FontSize, ...
              'FontName', 'cmr10', 'FontWeight', 'normal', 'FontAngle', 'normal' );
%% --------------------------------------------------------------------- %%
end