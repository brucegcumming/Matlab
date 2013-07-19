function result = PlotSims(data, varargin)
%PlotSims(data, varargin) plots results from dori simulations

if isfield(data,'dori') & length(data.dori) > 1
    plotmode = 1;
elseif isfield(data,'dsf')
plotmode = 2;
else
    plotmode = 1;
end


if plotmode == 1 % difference between +- dori, as a fucntino of dori, ori
ados = unique(abs(data.dori));
for k = 1:length(data.oris)
for j = 1:length(ados)
    id = find(abs(data.dori) == ados(j));
    if length(id) == 2
        sig(k,j) = diff(data.resps(k,id));
    end
end
end
plot(sig','o-');
result.sig = sig;
    result.doris = ados;
elseif plotmode == 2
asfs = unique(abs(data.dsf));
for k = 1:length(data.oris)
for j = 1:length(asfs)
    id = find(abs(data.dsf) == asfs(j));
    if length(id) == 2
        sig(k,j) = diff(data.resps(k,id));
    end
end
end
    result.dsfs = asfs;
plot(sig','o-');
result.sig = sig;
end