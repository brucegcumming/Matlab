function result = PlotRCLFPAll(rc, Expt, varargin)
%
%

%TO DO:
%subtract mean overall LFP to remove (a) framerate, (b) onset transient
%get peak time for spike var for the bystim analysis

trigbystim = 0; % set this to 1 to test Spike triggered lfp in viciniyt of stim
colors = mycolors;
lfptag = 'LFPRC';
j = 1;
submean = 0;
nopsych = 0;

while j <= length(varargin)
    if strncmpi(varargin{j},'bystim',4)
        trigbystim = 1;
    elseif strncmpi(varargin{j},'spk',3)
        result = PlotSpkMeans(rc);
        return;
    elseif strncmpi(varargin{j},'chan',3)
        j = j+1;
        lfpch = varargin{j};
    elseif strncmpi(varargin{j},'resp',3)
        j = j+1;
        stim = varargin{j};
        result = PlotResps(rc,stim);
        return;
    elseif strncmpi(varargin{j},'submean',4)
        submean = 1;
    elseif strncmpi(varargin{j},'tuning',3)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            trange = varargin{j};
        else
            trange = [900 1400]
        end
        result = PlotTuning(rc,trange);
        return;
    end
    j = j+1;
end
GetFigure(lfptag);


nch = size(Expt.Trials(1).LFP,2);
for j = 1:nch
    results{j} = PlotRCLFP(rc, Expt, 'chan',j,'submean');
    blanks(j,:) = results{j}.extralfp(1,:);
end

result.blanks = blanks;
result.results = results;


function avgs = PlotSpkMeans(res)

spk = res.results{1}.davg;
for j = 1:size(spk,2)
    mn(j) = mean(spk([1:10 end-10:end],j));
    avgs(:,j) = spk(:,j) - mn(j);
end
imagesc(avgs);

function resps = PlotTuning(res, times)

for j = 1:length(res.results)
%    id = find(res.results{j}.times >= times(1) & res.results{j} <= times(2);
    id = times(1):times(2);
    for k = 1:size(res.results{j}.lfpall,1)
        resps(j,k) = mean(res.results{j}.lfpall(k,id));
    end
end
imagesc(resps);

function resps = PlotResps(res, stim)

for j = 1:length(res.results)
%    id = find(res.results{j}.times >= times(1) & res.results{j} <= times(2);
    if stim == 0
        mn(j) = mean(res.results{j}.meansdf([1:10 end-10:end]));
        
        resps(j,:) = res.results{j}.meansdf - mn(j);
    else
        resps(j,:) = res.results{j}.lfpall(stim,:);
    end
end
imagesc(resps);

