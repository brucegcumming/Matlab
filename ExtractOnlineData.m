function Data = ExtractOnlineData(name, varargin)

[Data, Rs] = ReadDir(name);

for j = 1:length(Rs)
    t(j) = Rs{j}.stimlvl.times(1);
end
[ts, tid] = sort(t);
stimlvl = Rs{tid(1)}.stimlvl;
for j = 2:length(tid)
    stimlvl.times = cat(1,stimlvl.times,Rs{tid(j)}.stimlvl.times);
    stimlvl.level = cat(1,stimlvl.level,Rs{tid(j)}.stimlvl.level);
end
Data.(Rs{1}.stimch) = stimlvl;

function [DATA, RawData] = ReadDir(name, varargin)  %% Online Plots

d = dir(name);
reindex = 0;
args = {};
j = 1;
RawData = {};
nraw = 0;


while j <= nargin-2
        if strncmpi(varargin{j},'relist',3)
            reindex =1;
            args = {args{:} 'relist'};
        else
            args = {args{:} varargin{j}};
        end
    j = j+1;
end
    expnames = {};

        nexp = 1;
        SpkId = [];
        Spikes = [];
        Trialids = [];
        TrialTimes = [];
    lastn = nexp;
    DATA.badnames = [];
    Defaults.starttrial = 1;
    Defaults.fz = 100;
    for j = 1:length(d)
        if regexp(d(j).name,'Expt[0-9]*.mat') & ...
                d(j).bytes > 128 & .....
                isempty(strfind(d(j).name,'idx.mat')) & ...    %exclude the .idx file
                d(j).datenum < now-(1/(24 * 60 * 60)) & ... %at least  1 sec old
                isempty(strmatch(d(j).name,{expnames{:}}))%dont read if we already have
            if length(DATA.badnames)
                id = strmatch((d(j).name),DATA.badnames(badidx));
            else
                id = [];
            end
            if isempty(id) || d(j).datenum > DATA.badtimes(id(end))
            [trls, exps, All, Raw] = APlaySpkFile([name '/' d(j).name],'Defaults',Defaults,'online',args{:});
            if ~isempty(Raw)
                nraw = nraw+1;
                RawData{nraw} = Raw;
            end
            if isempty(exps)
                fprintf('%s No expts\n',d(j).name);
                DATA.badnames{nbad} = d(j).name;
                DATA.badtimes(nbad) = d(j).datenum;
                nbad = nbad+1;
                badidx = 1:nbad-1;
            else
                Expts{nexp} = exps{1};
                Expts{nexp}.gui.classified = 0;
                Expts{nexp}.gui.counted = 0;
                Expts{nexp}.gui.clustertype = 0;
                Expts{nexp}.gui.firsttrial = 1+length(Trialids);
                SpkId = [SpkId; trls.Spkid];
                newt = [Expts{nexp}.Trials.Trial];
                Trialids = [Trialids newt];
                for k = 1:length(Expts{nexp}.Trials)
                    news(k) = Expts{nexp}.Trials(k).Start(1);
                end
                TrialTimes = [TrialTimes news];
                Expts{nexp}.gui.ntrials = length(newt);
                if nexp > 1 && isfield(Spikes,'values') && size(Spikes.values,2) == size(All.Spikes.values,2)
                    Spikes.values = [Spikes.values; All.Spikes.values];
                    if isfield(All.Spikes,'dVdt')
                    Spikes.dVdt = [Spikes.dVdt; All.Spikes.dVdt];
                    end
                    Spikes.codes = [Spikes.codes; All.Spikes.codes];
                    Spikes.times = [Spikes.times; All.Spikes.times];
                else
                    Spikes = All.Spikes;
                end

                nexp = nexp+1;
                Defaults.starttrial = 1+ exps{1}.Trials(end).Trial;
                
                         if isfield(trls,'Probes');
                             for j = 1:length(trls.Probes)
                                 probenames{trls.Probes(j).probe} = trls.Probes(j).var;
                                 Probes(trls.Probes(j).probe) = trls.Probes(j);
                             end
                         end

            end
        end
        end
    end

    if exist('probenames','var') %ca be empty for online relist
    np = 0;
    probes = [];
    for j = 1:length(probenames)
        if ~isempty(probenames{j})
            np = np+1;
            DATA.probes(np) = Probes(j);
            probes = [probes j];
        end
    end
                DATA.probelist = probes;
            DATA.probevars = {probenames{probes}};
            DATA.probenames = cellstr(int2str(probes'));

    end


    if nexp == 1
        questdlg(sprintf('No expts in %s',name),'test','OK','OK');
    end
    DATA.name = name;
    DATA.state.online = 1;
    

