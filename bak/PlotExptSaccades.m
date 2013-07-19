function Saccades = PlotExptSaccades(Expt, varargin)

MARKTIMES = 1;
plottype = MARKTIMES;
minamp = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'minamp',6)
        j = j+1;
        minamp = varargin{j};
    end
    j = j+1;
end

ns = 1;

for j = 1:length(Expt.Trials)
   if length(Expt.Trials(j).Saccades)
   id = find([Expt.Trials(j).Saccades.size] > minamp);
   for k = 1:length(id)
       Saccades(ns) = Expt.Trials(j).Saccades(id(k));
       ns = ns+1;
   end
   end
end

if plottype == MARKTIMES
    ts = [[Saccades.peakt]; [Saccades.peakt]];
    ys = zeros(size(ts));
    yl = get(gca,'ylim');
    ys(2,:) = yl(1)+(diff(yl))/4;
    ys(1,:) = yl(1);
    line(ts, ys);
end