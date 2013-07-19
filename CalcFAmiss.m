function [famiss, details] = CalcFAmiss(Expt)

exa = Expt.Stimvals.et;
exb = Expt.Stimvals.e2;
id = find(abs([Expt.Trials.RespDir]) ==1);
Trials = Expt.Trials(id);
avals = [Trials.(exa)];
bvals = [Trials.(exb)];
choices =  [Trials.RespDir];
if isfield(Trials,'count')
    counts = [Trials.count];
elseif isfield(Trials,'Count')
    counts = [Trials.Count];
end
av = unique(avals);
bv = unique(bvals);
nmin = 4;

nx = 0;
for j = 1:length(av)
    for k = 1:length(bv)
        if abs(av(j)) > 0.0001
            if k == 1
                id = find(avals ==av(j));
                rates(j) = mean(counts(id));
                choicerate(j) = mean(choices(id));
            else
                id = [];
            end
        else
            id = find(avals ==av(j) & bvals == bv(k));
            id = []; %only use signal trials for now
        end
        if sum(choices(id) == -1) > nmin && sum(choices(id) == 1) > nmin;
            nx = nx+1;
            idxs{nx} = id;
            zscores{nx} = (counts(id)-mean(counts(id)))./std(counts(id));
            zchoices{nx} = (choices(id) ==1); %0 or 1
            details.zcp(nx) = CalcCP(counts(id(choices(id) ==1)),counts(id(choices(id) ==-1)));
            pp(nx) = sum(zchoices{nx})./length(zchoices{nx});
            stimval(nx) = av(j);
        end
        ns(j,k) = length(id);
    end
end

if nx < 2
    famiss = NaN;
    details.percentcorrect = 0;
    details.famiss = 0;
    return;
end
details.ti = (mean(rates(av > 0))-mean(rates(av < 0)))./(mean(rates(av < 0))+mean(rates(av > 0)));
details.pi = (mean(choicerate(av < 0))-mean(choicerate(av > 0)));
%pi > 0 means + stim values are associated with choice -1, so percentcorrect
%will be < 0.5. If ti is -v, this will be correct = p (detecting) rate
%increase is lower 
for j = 2:nx
    for k = 1:j-1 
        pc(j,k) = (sum(zchoices{k} ==0) + sum(zchoices{j} ==1))./(length(zchoices{j}) + length(zchoices{k}));
        nrates = (counts(idxs{k}) - mean(counts(idxs{k})))./(mean(counts(idxs{j}))-mean(counts(idxs{k})));
        prates = (counts(idxs{j}) - mean(counts(idxs{k})))./(mean(counts(idxs{j}))-mean(counts(idxs{k})));
        famiss(j,k) = mean(nrates(zchoices{k} == 1))-mean(prates(zchoices{j} ==0));
    end
end
id = find(pc > 0);
details.percentcorrect = pc(id);
details.famiss = famiss(id);
famiss = mean(famiss);
