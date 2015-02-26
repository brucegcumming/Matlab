function CombineAll(a,DATA, varargin)
%
%Combine all probes, for one expt
%Don't redo LFP here. Too slow.

checkonce = 1;
guicall = 0;
cellsonly = 0;

j = 1;

dolfp = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'check',3)
        checkonce = 1;
    elseif strncmpi(varargin{j},'lfp',3)
        dolfp = 1;
    elseif sum(strncmpi(varargin{j},{'cells' 'mucells'}, 5))
        cellsonly = 1;
        cmb.CombineAllCells(a, DATA, varargin{:});
        return;
    elseif strncmpi(varargin{j},'nocheck',3)
        checkonce = 0;
    end
    j = j+1;
end

if ~isstruct(DATA)
    DATA = GetDataFromFig(a);
    guicall = 1;
end
DATA.state.includeprobename = 1;
eid = get(DATA.clst,'value');
outname = DATA.outname;
DATA.extype = eid;
oldauto = DATA.plot.autoclustermode;
%hitting "combine all" online is asking to autocut uncut probes
if DATA.plot.autoclustermode == 3 && DATA.state.online
    DATA.plot.autoclustermode = 1;
end

%Make sure all the spike times are loaded first.
%this makes the spike loading quicker when online
if DATA.state.online
    ids = get(DATA.elst,'value');
   DATA =  cmb.BuildAllTimes(DATA,'expts',DATA.expid(ids));
   SetData(DATA);
end
ts = now;
sumload = 0;
oldplot = DATA.plot;
oldstate = DATA.state;
DATA.plot.quickspks = 0;
DATA.state.showspikes = 0;
for j = 1:length(DATA.probelist)
    DATA.loaddur = 0;
    DATA = cmb.SetProbe(DATA, DATA.probelist(j));
    loaddurs(j,:) = DATA.loaddur;
    if j == checkonce
        chk = DATA.state.interactive;
    else
        chk = 0;
    end
    [Expt, a, plotres] = cmb.CombinePlot(DATA, chk);
%DATA.spkcache could have been changed so copy this.
%This is a place where a separate appdata would work better
    DATA = CopyFields(DATA, a,'spkcache');
    Expt.plotres = rmfield(plotres,'Data');
    Expt.Header.probe = DATA.probe;
    if ~isfield(Expt.Header,'cellnumber') && DATA.state.online 
        Expt.Header.cellnumber = 0;
    end
    if isfield(DATA,'AllClusters') && ~isfield(Expt.Header,'probesep')
        Expt.Header.probesep = 50;
    end
    %    Expt.Header.SpkStats = cmb.GetSpkStats(DATA);
    drawnow;
    nspk = sum([Expt.Trials.count]);
    file = regexprep(outname,'\.p[0-9]*c1\.',['.p' num2str(DATA.probelist(j)) 'c1.']);
    %    file = cmb.CombinedName(DATA,eid,1);
    Expts{j} = Expt;
    if (DATA.state.nospikes == 0 || DATA.state.nospikes == 2) && DATA.state.online == 0
        BackupFile(file,'print');
        Expt = rmfields(Expt,'plotres');
        save(file,'Expt');
        fprintf('Saved %d spikes to %s\n',nspk,file);
    end
    
    if isfield(Expt,'suffixes')
        estr = sprintf('%s Suffs %s',sprintf(' %d',Expt.Header.Combineids),sprintf(' %d',Expt.Header.suffixes));
    else
        estr = sprintf(' %d',Expt.Header.Combineids);
    end
    if DATA.logfid
        fprintf(DATA.logfid, '%s,Saved %d spikes (Expts%s) to %s (%s) (all%d)\n',datestr(now),nspk,estr,file,DATA.user,j);
    end
end
DATA.Expt = Expt;
DATA.plot.autoclustermode = oldauto;
fprintf('Building All Expts took %.2f (%.2f,%.2f loading spikes)\n',mytoc(ts),sum(loaddurs));

if DATA.state.online == 0 && dolfp
    cmb.combine('savelfp',DATA);
end
DATA.AllExpts = Expts;
DATA.plot.autoclustermode = oldauto;
DATA.plot = CopyFields(DATA.plot,oldplot,'quickpsks', 'autoclustermode');
DATA.state = CopyFields(DATA.state,oldstate,'showspikes');
if DATA.state.interactive || guicall
    set(DATA.toplevel,'UserData',DATA);
    PlotAllCellFiles(DATA.AllExpts,'parentfigure',DATA.toplevel,'tag',DATA.tag.allexpts);
elseif guicall
    set(DATA.toplevel,'UserData',DATA);
end

