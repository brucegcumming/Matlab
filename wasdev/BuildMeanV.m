function sumv = BuildMeanV(datadir, e, probes, varargin)
%Build MeanV file from individual Grid probe files.
%sumv = BuildMeanV(datadir, expt, probes, varargin)
%...,'save')  Saves result.
% Called by BuildGridFullV.
savefile = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'save',4)
        savefile = 1;
    end
    j = j+1;
end
    for j = 1:length(probes)
        p = probes(j);
        load(sprintf('%s/Expt%d.p%dFullV.mat',datadir,e,p));
        V = double(FullV.V);
        if j == 1
            sumv = V./std(V);
        else
            sumv(1:length(V)) = sumv(1:length(V)) + V./std(V);
        end
    end
    sumv = sumv./length(probes);
    MeanV.savetime = now;
    if savefile
        save(sprintf('%s/Expt%dFullVmean.mat',datadir,e),'sumv','MeanV')
    end
