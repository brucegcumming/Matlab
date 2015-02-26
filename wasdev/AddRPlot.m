function [rax,h] = AddRPlot(lax, varargin)
%[rax,h] = AddRPlot(lax, varargin)
% add an additional axis on the RHS of a graph which already has an axis (lax) on
% the left

if isfigure(lax)
    fig = lax;
    axs = get(fig,'Children');
    lax = axs(1);
else
fig = get(lax,'Parent');
axs = get(fig,'Children');
end


set(lax,'Box','off');
ax1hv = get(lax,'HandleVisibility');
rax = [];
for j = 1:length(axs)
    p = get(axs(j),'YaxisLocation');
    if strcmp(p,'Right')
        rax = axs(j);
    end
end

if isempty(rax) || force
rax = axes('HandleVisibility',ax1hv,'Units',get(lax,'Units'), ...
    'Position',get(lax,'Position'),'Parent',fig);
end
%%raxis = axes('Position',get(gca,'Position'),'color','None');
if nargin > 1
    h = plot(varargin{:});
end
set(rax,'YAxisLocation','right','Color','none', .... 
          'XGrid','off','YGrid','off','Box','off','XtickLabel',[],...
          'Xlim',get(lax,'Xlim'),'Xscale',get(lax,'Xscale'));
