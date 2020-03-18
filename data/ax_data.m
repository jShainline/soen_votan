function [x,y] = ax_data(ax,n)
% gets 'x' and 'y' data for curve number 'n' from axis 'ax'
h = get(ax,'children');
x = get(h(length(h)+1-n),'xdata');
y = get(h(length(h)+1-n),'ydata');
