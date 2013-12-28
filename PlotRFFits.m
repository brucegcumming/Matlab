function PlotRFFits(fits, varargin)
%plots RF fits made with BuildRFFits
colors = mycolors;
np=0;
hold off;
for j = 1:length(fits)
    if isfield(fits{j},'proberf') && size(fits{j}.proberf,1) > 2
        np = np+1;
        depths = fits{j}.proberf(:,11)+fits{j}.proberf(:,12) .* fits{j}.spacing./1000;
        [depths,b] = sort(depths);
        rfs = (fits{j}.proberf(:,1)-fits{j}.proberf(:,9)) +i*(fits{j}.proberf(:,2)-fits{j}.proberf(:,10));
        nc = 1+mod(np-1,10);
        plot(rfs(b),'o-','color',colors{nc},'buttondownfcn',{@HitPen, j});
        hold on;
    end
end
setappdata(gcf,'fits',fits)

function HitPen(a,b, pen);
F = GetFigure(a);
fits = getappdata(F,'fits');
fprintf('Fit %d:%s\n',pen,fits{pen}.dirname);
j = pen;
depths = fits{j}.proberf(:,11)+fits{j}.proberf(:,12) .* fits{j}.spacing./1000;
[depths,b] = sort(depths);
rfs = (fits{j}.proberf(:,1)-fits{j}.proberf(:,9)) +i*(fits{j}.proberf(:,2)-fits{j}.proberf(:,10));
rfs = rfs(b);
for j = 1:length(rfs)
    fprintf('%.3f %.2f %.2f\n',depths(j),real(rfs(j)),imag(rfs(j)));
end
   
