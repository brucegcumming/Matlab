function Expt = SetPdir(Expt,varargin)
verbose = 0;
bysign = 0;

j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'bysign',5)
    bysign = 1;
elseif strncmpi(varargin{j},'verbose',4)
    verbose = 1;
end
j = j+1;
end
if isfield(Expt.Trials,'dfx')
    dfx = [Expt.Trials.dfx] - [Expt.Trials.fx];
else
    dfx = zeros(1,length(Expt.Trials));
end
if isfield(Expt.Trials,'dfy')
    dfy = [Expt.Trials.dfy] - [Expt.Trials.fy];
else
    dfy = zeros(1,length(Expt.Trials));
end
Fa = (dfx + i * dfy);
or = GetEval(Expt,'or');
d = cos(angle(Fa) + pi/2 - or .*pi/180) .* abs(Fa);
d = 10000 .* d./[Expt.Trials.dur];
d = round(d .* 5)./5;
for j = 1:length(Expt.Trials)
    if bysign
        Expt.Trials(j).pvel = sign(d);
    else
        Expt.Trials(j).pvel = d;
    end
end

Expt.Stimvals.Pursuitdir = 90 + mean(angle(Fa(find(sign(d) > 0)))) * 180/pi;

if verbose
    fprintf('Ori %.0f\n',GetEval(Expt,'or'));
ds = unique(d);
for j = 1:length(ds)
    id = find(d == ds(j));
    df(j) = mean(Fa(id));
    fprintf('%.2f: %.2f %.2f\n',ds(j),real(df(j)),imag(df(j)));
end
end