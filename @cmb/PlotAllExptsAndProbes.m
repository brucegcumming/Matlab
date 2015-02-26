function PlotAllExptsAndProbes(a,b, varargin)
if isfield(a, 'state') && isfield(a,'AllData')
DATA = a;
else
DATA = GetDataFromFig(a);
end
tstart = now;

DATA.state.autoplotnewprobe = 0;

for j = 1:length(DATA.probelist)
DATA = cmb.SetProbe(DATA, DATA.probelist(j));
cmb.PlotAllExpts(DATA, 0, 'figlabel', ['AllExpts Probe ' num2str(j)]);
end
mytoc(tstart);
DATA.state.autoplotnewprobe = 1;

