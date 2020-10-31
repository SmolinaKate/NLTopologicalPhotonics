function [ position ] = ExpandAxesToMinimizeWhiteSpace( ax )
%% --------------------------------------------------------------------- %%
%% Expand Axes to Fill Maximum Available Space and Minimize White Space.
% The input parameter is an axes object.
% The output parameter is a four-element vector with size and location an axes object.
%% --------------------------------------------------------------------- %%
ax.Units = 'normalized';
ax.OuterPosition = [ 0, 0, 1, 1 ];
%% --------------------------------------------------------------------- %%
% Get the dimensions of the maximum available space from the OuterPosition property of the axes.
outer_position = ax.OuterPosition;
% Account for the space needed for the tick values and text labels using the margin values stored
% in the TightInset property (margins for the text labels). Notice this property is read-only.
tight_inset = ax.TightInset;
%% --------------------------------------------------------------------- %%
left = outer_position(1) + tight_inset(1);
bottom = outer_position(2) + tight_inset(2);
width = outer_position(3) - tight_inset(1) - tight_inset(3);
height = outer_position(4) - tight_inset(2) - tight_inset(4);
position = [ left, bottom, width, height ];
%% --------------------------------------------------------------------- %%
% Expand the axes size so that it fills the maximum available space in the figure.
ax.Position = position;
%% --------------------------------------------------------------------- %%
end