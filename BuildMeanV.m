function varargout = BuildMeanV(datadir, e, probes, varargin)
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

if nargin == 1
    BuildMeanForDir(datadir);
    return;
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
    MeanV.probes = probes;
    if savefile
        save(sprintf('%s/Expt%dFullVmean.mat',datadir,e),'sumv','MeanV')
    end
    varargout{1} = sumv;
    
 function    BuildMeanForDir(datadir)

     
     d = dir([datadir '/*FullV.mat']);
     for j = 1:length(d)
         expts(j) = GetExptNumber(d(j).name);
         probes(j) = GetProbeFromName(d(j).name);
     end
     Array = GetArrayConfig(datadir);
     probes = unique(probes);
     if isfield(Array,'badprobes')
         probes = setdiff(probes,find(Array.badprobes));
     end
     expts = unique(expts);
     for j =1:length(expts)
         fullvname = [datadir '/Expt' num2str(expts(j)) 'FullVmean.mat'];
         if ~exist(fullvname)
             fprintf('Making %s\n',fullvname);
             BuildMeanV(datadir, expts(j), probes, 'save');
         else
             fprintf('%s already exists\n',fullvname);
         end
     end