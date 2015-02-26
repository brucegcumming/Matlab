function DATA = AutoCut(DATA, eid, idx, varargin)
%DATA = AutoCut(DATA, eid, varargin)
plotxy =0;
recalc = 0;
automode = 'energy';
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'recalc')
       recalc = 1;
    elseif strcmp(varargin{j},'show')
        plotxy = 1;
    end
    j = j+1;
end

if DATA.state.online
%    DATA = cmb.LoadSpikes(DATA,eid);
end
%DATA = cmb.AutoCut(DATA,eid,idx,'cluster',varargin{:});
%use previous autocut if possible
e = eid;
p = DATA.probe;

if recalc == 0 && cmb.SpkCache(DATA, e, p, 'check')
    X = cmb.SpkCache(DATA, e, p, 'get');
    idx = find(X.codes ==1);
    fprintf('Using Existing codes forE%dP%d\n',e,p);
    automode = 'none';
elseif DATA.state.online && DATA.AllData.Spikes.exptid ~= e
    DATA = cmb.LoadSpikes(DATA,eid);
end
ts = now;
if isempty(DATA.spklist) && ~strcmp(automode,'none')
    DATA = cmb.SpkCache(DATA,e,p,'set');
    return;
end

if strcmp(automode,'pcs')
    fprintf('Using PCA to autocut E%dP%d',e,p);
    [V,E] = eig(cov(DATA.AllData.Spikes.values));
    scores = DATA.AllData.Spikes.values * V;
    if E(end) > E(1)
        scores = fliplr(scores);
    end
    D = MyDip(scores(:,1));
    crit = D.x(D.dip(1));
    idx = find(scores(:,1) > crit);
    fprintf('   took %.2f\n',mytoc(ts));
elseif strcmp(automode,'energy')
    fprintf('Using Energy to autocut E%dP%d',e,p);
    idx = find(DATA.Spikes.cx(DATA.spklist) ~= 0);
    D = MyDip(DATA.Spikes.cx(DATA.spklist(idx)),'quick');
    dips.npeaks = length(D.peaks);
    dips.sm = D.sm(end);
    if length(D.sm) > 1
        dips.shift = sign(diff(D.sm(1:2)));
    else
        dips.shift = 0;
    end
    dips = CopyFields(dips,D,'sigma','nvals','peakcounts');
    DATA.quickdips(e,p) = dips;
    if ~isempty(DATA.spklist)
    crit = D.x(D.dip(1));
    idx = find(DATA.Spikes.cx(DATA.spklist) > crit);
    a = abs(crit);
    while isempty(idx) && abs(crit) > a/4
        crit = crit * 0.9;
        idx = find(DATA.Spikes.cx(DATA.spklist) > crit);
    end
    
    idx = DATA.spklist(idx);
    end
    fprintf('   took %.2f\n',mytoc(ts));
end
DATA.AllData.Spikes.codes(:,2) = 0;
DATA.AllData.Spikes.codes(idx,2) = 1;
if ~strcmp(automode,'none')
    DATA = cmb.SpkCache(DATA, eid, DATA.probe,'set');
end
if plotxy
    cmb.SetFigure(DATA.tag.clusterxy,DATA);
    hold off;
    cmb.DrawXYPlot(DATA, DATA.spklist);
end