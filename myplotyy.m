function myplotyy(x,ya, yb, varargin)

     plot(x,ya);
raxis = axes('Position',get(gca,'Position'));
plot(x,yb);
set(raxis,'Color','None','Xlim',get(gca,'Xlim'),'YAxisLocation','Right','Xscale','Log');

