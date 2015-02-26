function PlotRFFits(fits, varargin)
%plots RF fits made with BuildRFFits
colors = mycolors;
plotbypen = 0;
np=0;
name = [];
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'bypen',4)
        plotbypen = 1;
    elseif strncmp(varargin{j},'name',4)
        j = j+1;
        name = varargin{j};
    end
    j=j+1;
end

if ~isempty(name)
    for j = 1:length(fits)
        if strncmp(name,GetName(fits{j}),length(name))
            good(j) = 1;
        else
            good(j) = 0;
        end
    end
    fits = fits(find(good));
end


if plotbypen
pe = CellToMat(fits,'Pn');
pes = unique(pe);
pes = pes(pes > 0);
for j = 1:length(pes)
    id = find(pe == pes(j));
    penfit{j} = fits{id(1)};;
    for k = 2:length(id)
        penfit{j}.proberf = cat(1,penfit{j}.proberf,fits{id(k)}.proberf);
    end
end
fits = penfit;
end
GetFigure('RF Fits');
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

GetFigure('Fit Data')
%PlotExptFit(fits{j}.fits);
   
