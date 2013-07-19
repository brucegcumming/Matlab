function MakeFigure(varargin)

resps = [];
tags{1} = 'DTsep';
pix2deg = 0.05;
fsz = 12; 

if length(varargin) > 1 & isstruct(varargin{1})
    details = varargin{1};
    resps = varargin{2};
    j = 3;
else
    j = 1;
end
stimr = 3;

[resps, details] = AbsSlant('rfp',[-1 -1],'dtscale',[0.4 3],'dxs',[-0.2 0 0.1 0.3 0.4],'even','dtrange',[-1:0.05:1]); 
iw = details.iw;
ih = details.ih;
GetFigure(tags{1});
subplot(3,2,3);
hold off;
plot(details.dtx,details.dtf,'k','linewidth',2);
colors = mycolors('matlab');
colors{5} = [0  0.75 0.75];
hold on;
for j = 1:length(details.disps)
plot(details.disps(j),max(details.dtf).*1.05,'v','color',colors{j},'markerfacecolor',colors{j});
end
hold off;
subplot(3,2,[1 2]);
hold off;
imagesc(details.rf);
axis('image');
set(gca,'XTicklabel',[], 'YTicklabel',[])

h = title('Receptive field and Stimulus');
set(h,'fontsiz',fsz);
hold on;
a = 0:pi./20:2*pi;
r = ones(size(a)) .* stimr./pix2deg;
[x,y] =pol2cart(a,r);
plot(x+iw/2,y+ih/2,'w--');
subplot(3,2,4);
hold off;
for j = 1:length(details.disps)
plot(details.slants .* 180/pi,resps(j,:),'color',colors{j},'linewidth',2);
hold on;
end
set(gca,'Xtick',[0 90 180 270 360],'xlim',[0 360]);

[resps, details] = AbsSlant('rfp',[-1 -1],'dtscale',[0.4 60],'dxs',[0 0.25 0.3 0.4 0.45],'dtrange',[-1:0.05:1]);
subplot(3,2,6);
hold off;
for j = 1:length(details.disps)
plot(details.slants .* 180/pi,resps(j,:),'color',colors{j},'linewidth',2);
hold on;
end
xlabel('Tilt Angle (degrees)','fontsiz',fsz);
set(gca,'Xtick',[0 90 180 270 360],'xlim',[0 360]);
subplot(3,2,5);
hold off;
plot(details.dtx,details.dtf,'k','linewidth',2);
hold on;
for j = 1:length(details.disps)
plot(details.disps(j),max(details.dtf).*1.05,'v','color',colors{j},'markerfacecolor',colors{j});
end
xlabel('Disparity','fontsiz',fsz);
hold off;

