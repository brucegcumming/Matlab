function spk = spktimes(E, varargin)
%return list of absolute spike times, in sec

if iscell(E)
    for j = 1:length(E)
        spk{j} = expt.spktimes(E{j},varargin{:});
    end
    return;    
elseif isfield(E,'Spikes') && isfield(E,'trialid') %Allexpt Spikes strcut
    spk = cat(1,E.Spikes{:})./10000;
    return;
elseif isfield(E,'Spikes') && iscell(E.Spikes) %AllExpt Not working yet
    for j = 1:length(E.Spikes)
        spks = [];
        for t = 1:length(E.Spikes{j}.Spikes)
            spks = [spks E.ExptTrials(t).Start(1)+ E.Spikes{j}.Spikes{t}'];
        end
        spk{j} = spks;
    end
    return;
else
end

spk = [];
    for j = 1:length(E.Trials)
        spk = [spk E.Trials(j).Start(1)+E.Trials(j).Spikes'];
    end
spk = spk ./10000;