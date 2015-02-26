function tid = FindTrial(Expt, t, varargin)
%tid = FindTrial(Expt, t, varargin) find Trial that contains time t

for j = 1:length(Expt.Trials)
    starts(j) = Expt.Trials(j).Start(1)./10000;
    ends(j) = Expt.Trials(j).End(end)./10000;
end
tid = find(t > starts & t < ends);
if isempty(t)
    [a,b] = min(abs(t-starts));
    [c,d] = min(abs(t-ends));
    if a > c
        tid = d;
    else
        tid = b;
    end
end
