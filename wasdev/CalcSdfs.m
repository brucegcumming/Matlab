function sdf = CalcSdfs(Expt, type, width)

latency = 500;
duration = min([Expt.Trials.End] - [Expt.Trials.Start]);
times = -260:20:duration+latency;
vallist = eval(['[Expt.Trials.' type ']']);
vals = unique(vallist);

[Expt.Trials.Trigger] = deal(0);
j = 1;
sdf.times = times;
for val = vals;
  idx = find(vallist == val);
  [tmp, n] = trigsdf(Expt.Trials(idx), width, times, 'exp');
  sdf.sdf(j,1:length(tmp)) = deal(tmp);
  sdf.val(j) = val;
  sdf.n(j) = n;
  j = j+1;
end

