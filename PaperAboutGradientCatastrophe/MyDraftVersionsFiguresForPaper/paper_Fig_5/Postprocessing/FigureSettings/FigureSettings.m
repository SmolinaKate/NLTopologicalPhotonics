function [ ] = FigureSettings( fig, papersize )
%% --------------------------------------------------------------------- %%
%% Position.
fig.Units = 'centimeters';                  % 'normalized' | 'inches' | 'points' | 'pixels'
fig.Position = [ papersize./ 2, papersize ];
fig.PaperUnits = 'centimeters';
fig.PaperSize = papersize;
fig.PaperPosition = [ 0, 0, papersize ];
%% --------------------------------------------------------------------- %%
fig.Color = [ 1, 1, 1 ];
%% --------------------------------------------------------------------- %%
end
