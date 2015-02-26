function Expt = AddExpts(varargin)
%Combines a list of Expts into one,
%Expt = AddExpts(expt1, expt2, expt3....)
%fills fields of Expt.Trials to be union of those in expt1 to exptn
%See also CombineExpts

Expt = [];

explist = varargin;

for j = 1:length(explist)
    slvals(j) = explist{j}.Stimvals.sl;
end

if length(unique(slvals)) > 1
    setseedlen = 1;
    seeds = explist{1}.Stimvals.sl .* ones(size(explist{1}.Trials(1).Start));
    [explist{1}.Trials.sl] = deal(seeds);
else
    setseedlen = 0;
end
    
Expt.Header = explist{1}.Header;
Expt.Stimvals = explist{1}.Stimvals;
Expt.Trials = explist{1}.Trials;
Expt.Names = explist{1}.Names;
Expt.isrc = explist{1}.isrc;

fnames = fieldnames(Expt.Trials);
for j = 1:length(Expt.Trials)
    slen(j) = length(Expt.Trials(j).Start);
end
startlen = ceil(prctile(slen,20));

nt = length(Expt.Trials)+1;
for j = 2:length(explist)
    if setseedlen
        seeds = explist{j}.Stimvals.sl .* ones(size(explist{j}.Trials(1).Start));
        [explist{j}.Trials.sl] = deal(seeds);
    end
    for trial = 1:length(explist{j}.Trials)
        if length(explist{j}.Trials(trial).Start) >= startlen
            for k = 1:length(fnames)
                if isfield(explist{j}.Trials,fnames{k})
                    Expt.Trials(nt).(fnames{k}) = explist{j}.Trials(trial).(fnames{k});
                else
                    Expt.Trials(nt).(fnames{k}) = explist{j}.Stimvals.(fnames{k});
                end
            end

            nt = nt+1;
        end
    end
%    Expt.Trials = [Expt.Trials explist{j}.Trials];
end

if length(unique(slvals)) > 1
    if strcmp(Expt.Stimvals.e2,'e0')
        Expt.Stimvals.e2 = 'sl';
    end
end
