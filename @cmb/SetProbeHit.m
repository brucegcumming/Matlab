function SetProbeHit(a,b, varargin)

next = 0;
DATA = GetDataFromFig(a);
DATA.state.savedvdt = 0;
if DATA.savedclusters == 0
    msgbox('YOu have not saved cluster parameters for last probe yet');
    return;
end
if DATA.playingspk
    set(findobj(DATA.svfig,'Tag','StopSpool'),'value',1);
    return;
end

j = 1;
while j <= length(varargin)
    if strcmp('ReloadProbe',varargin{j})
        if isfield(DATA,'AllSpikes')
            DATA = rmfield(DATA,'AllSpikes');
        end
        if isfield(DATA,'AllClusters')
            DATA = rmfield(DATA,'AllClusters');
        end
        if isfield(DATA,'sids')
            DATA = rmfield(DATA,'sids');
        end
    elseif strcmp(varargin{j},'next')
        next = 1;
    elseif strcmp(varargin{j},'prev')
        next =  -1;
    end
    j = j+1;
end

f = get(a);
if strcmp(f.Tag,'ProbeId')
    pit = a;
else
    pit = findobj(DATA.toplevel,'Tag','ProbeId');
    sit = findobj(DATA.toplevel,'Tag','SubprobeId');
end
DATA.oldprobe = DATA.probe;
if isfield(f,'Tag') & strcmp(f.Tag,'SubprobeId')
    DATA.subprobe = f.Value-1;
elseif isfield(f,'Value')
    id = get(pit,'value');
    probelist = getappdata(pit,'probelist');
    if isempty(probelist)
        probelist = DATA.probelist;
    end
    if next
        newid = id+next;
        if (DATA.subprobe == 0 || isempty(sit)) && newid <= length(probelist) && newid > 0
            id = newid;
            set(pit,'value',id);
        elseif length(sit) && DATA.subprobe > 0
            sub = get(sit,'value');
            if sub < 5
                sub = sub+1;
                set(sit,'value',sub);
                DATA.subprobe = sub-1;
            elseif id < length(probelist)
                id = id+1;
                DATA.subprobe = 1;
                set(sit,'value',DATA.subprobe+1);
                set(pit,'value',id);
            end
        end
    end
    DATA.probe = probelist(id);
end


DATA = cmb.CheckState(DATA);

if isfield(DATA,'AllSpikes')
    if DATA.state.applylastcluster
        [DATA, DATA.spklist] = SetExptSpikes(DATA,DATA.currentexpt(1),0,'useexpt');
    else %?? Should't need cache with AllSpikes
        DATA = cmb.SpkCache(DATA, DATA.currentexpt{1},DATA.probe,'add');
        [DATA, DATA.spklist] = SetExptSpikes(DATA,DATA.currentexpt(1),'setrange');
    end
    nloaded = 0;
    for j = 1:length(DATA.AllSpikes)
        if isfield(DATA.AllSpikes{j},'codes') & length(DATA.AllSpikes{j}.codes) > 10
            nloaded = nloaded+1;
        end
    end
elseif DATA.probe == 100
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif isfield(DATA,'AllClusters')
    DATA.state.nospikes =2;
    DATA = cmb.SetProbe(DATA, DATA.probe);
else
    DATA = cmb.SetProbe(DATA, DATA.probe);
    nloaded = 1;
end
if DATA.listbycell
    DATA.outname = regexprep(DATA.outname,'.p[0-9]*c([0-9]).','.cell0.');
    DATA.outname = regexprep(DATA.outname,'.cell[0-9]*.',sprintf('.cell%d.',DATA.probe));
else
    DATA.outname = regexprep(DATA.outname,'.cell[0-9]*.','.p0c1.');
    DATA.outname = regexprep(DATA.outname,'.p[0-9]*c([0-9]).',sprintf('.p%dc$1.',DATA.probe));
end
set(DATA.saveitem,'string',DATA.outname);
set(DATA.toplevel,'UserData',DATA);
playspk = cmb.GetCheck('ShowSpikes');
DATA.state.showspikes = playspk;
playspk = (playspk | DATA.plot.quickspks);
%
% if > 3 probes are laoded, don't spool through them every time the probe
% hit is changed.
if DATA.state.nospikes
    fprintf('Calling combine from Setprobehit\n');
    DATA = cmb.combine('setexp', DATA,'newprobe');
    if DATA.state.verbose fprintf('Returned from  combine in Setprobehit\n'); end
    cmb.plotISI(DATA);
elseif (DATA.state.autospool | playspk) & nloaded < 4
    DATA = cmb.combine('setexp', DATA,'newprobe');
elseif DATA.state.usexycache
    DATA = cmb.combine('setexp', DATA,'newprobe');
else
    DATA = CalcClusterVars(DATA,  DATA.spklist);
    DATA.Expts{DATA.currentexpt(1)}.gui.spks = DATA.spklist;
    if isfigure(DATA.xyfig)
        GetFigure(DATA.xyfig);
        hold off;
        cmb.DrawXYPlot(DATA, DATA.spklist);
    end
end
set(DATA.toplevel,'UserData',DATA);


