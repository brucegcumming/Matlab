function examples(type, varargin)

if strncmpi(type,'plotyy',5)
    testplot;
end



function testplot

%figure;
fig = gcf;
x = 1:100;
y = rand(size(x));
plot(x,y);
hold on;
laxis = gca;
set(laxis,'color','none');
ax1hv = get(laxis,'HandleVisibility');
raxis = axes('HandleVisibility',ax1hv,'Units',get(laxis,'Units'), ...
    'Position',get(laxis,'Position'),'Parent',fig);
%%raxis = axes('Position',get(gca,'Position'),'color','None');
plot(x,y./10,'r');
set(raxis,'YAxisLocation','right','Color','none', ...
          'XGrid','off','YGrid','off','Box','off');

%set(raxis,'Color','None','Xlim',get(laxis,'Xlim'),'YAxisLocation','Right','Xscale','log');
%drawnow;
%set(raxis,'Xscale','lin');
%set(raxis,'color','none');

