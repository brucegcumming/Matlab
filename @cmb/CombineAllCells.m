function AllExpt = CombineAllCells(a,DATA, varargin)
dolfp = 0;
savefiles = 0; %save individual files.  Always makes AllExpt
makemu = 0;
guicall = 0;
oneprobe = 0;
useguilist = 1;

if ~isstruct(DATA)
    DATA = GetDataFromFig(a);
    guicall = 1;
end

j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'saveonlyall')
        savefiles = 0;
    elseif strcmp(varargin{j},'save')
        savefiles = 1;
    elseif strcmp(varargin{j},'mucells')
        makemu = 1;
    elseif strcmp(varargin{j},'muonly')
        makemu = 2;
    elseif strcmp(varargin{j},'oneprobe')
        j = j+1;
        oneprobe = varargin{j};
    elseif strcmp(varargin{j},'quick')
        DATA.state.interactive = 0; %don't plot graphs
    elseif strcmp(varargin{j},'reapply')
        useguilist = 0;
    end
    j = j+1;
end
if useguilist
    combinelist = get(DATA.elst,'value');
end

ts = now;
DATA.state.recount = 1;

addedid = [];

DATA.state.includeprobename = 1;
oldlistbycell = DATA.listbycell;

%should have a list of files by here. But if not will need to read list
%from each saved file.
useeachcombineids = 0;
checkonce =1;
eid = get(DATA.clst,'value');
outname = get(DATA.saveitem,'string');
%if the user has changed the filename by adding chars before the expt label
%(e.g.  grating.sf1OXM instead of grating.OXM), then find the modifier part
%and keep it.
modifier = '';
DATA.extype = eid;
file = cmb.CombinedName(DATA,DATA.extype,1);

id = regexp(file,'\.');
ndot = 2; %Needs to be 2 e.g. lemM214/ ? 3 at othe rtimes?
if length(id) > ndot
    expname =file(id(ndot)+1:id(ndot+1)-1);
end
id = regexp(outname,'\.');
if length(id) > ndot
    bexpname = outname(id(ndot)+1:id(ndot+1)-1);
    id = strfind(bexpname,expname);
    if ~strcmp(bexpname,expname) && ~isempty(id)
        modifier = bexpname(1:id-1);
    end
end

allexptname = regexprep(outname,'\.p[0-9]*c1\.','.Cells.');
allexptname = regexprep(allexptname,'\.c[0-9]*\.','.Cells.');
allexptname = regexprep(allexptname,'\.cell[0-9]*\.','.Cells.');


oldplot = DATA.plot;
oldauto = DATA.plot.autoclustermode;
oldflip = DATA.plot.flip;
oldspikes = DATA.state.showspikes;
DATA.state.showspikes = 0;  %don't need to load spike waveforms.
%hitting "combine all" online is asking to autocut uncut probes

%if saving files, do probe 1 first
DATA.listbycell = 0;
DATA.currenexpt = DATA.expid(1);
DATA = cmb.SetProbe(DATA, 1);
file = cmb.CombinedName(DATA,DATA.extype,1,'modifier',modifier);
if useguilist == 0
    H = [];
    if exist(allexptname)
        X = matfile(allexptname);
        if sum(strcmp('ExptHeader',fields(X)))
            H = X.ExptHeader;
        end
    end
    if isempty(H) && exist(file)
        load(file);
        H = Expt.Header;
    end
    if ~isempty(H)
%use exptnos not expid here - previous file may have manually combined over
%expt types
        combinelist = find(ismember(DATA.exptnos,H.Combineids));
        DATA.expid = combinelist;
        combinelist = 1:length(combinelist);
        [Expt, DATA, plotres] = cmb.CombinePlot(DATA, DATA.state.interactive,'ids',1:length(combinelist));
        Expt.Header.expid = combinelist;
    else
        combinelist = [];
        useeachcombineids = 1;
        Expt = [];
    end
else
    DATA.newload = 0;
    [Expt, DATA, plotres] = cmb.CombinePlot(DATA, DATA.state.interactive);
    if DATA.newload
        DATA.newload = 0;
        SetData(DATA);
    end
end
if isempty(Expt)
    return;
end
ExptHeader = Expt.Header;
nspk = sum([Expt.Trials.count]);
if isempty(combinelist)
    useeachcombineids = 1;
else %alwasy save out p1 because that records the combineids...
    ts = now;
    BackupFile(file,'print');
    save(file,'Expt','ExptHeader');
    fprintf('Saved %d spikes (Expts %s) to %s (%s) took %.2f\n',nspk,sprintf(' %d',DATA.combineids),file,DATA.user,mytoc(ts));
end
AllExpt.Expt = Expt;

if useguilist == 0
    ids = combinelist;
else
    ids = get(DATA.elst,'value');
end
if isfield(DATA,'CellDetails') && isfield(DATA.CellDetails,'exptids')
    eid = find(ismember(floor(DATA.CellDetails.exptids),DATA.exptnos(DATA.expid(ids))));
    cells = unique(DATA.CellList(eid,:,:));
    cells = cells(cells > 0);
    probelist = 1:size(DATA.CellList,2);
elseif isfield(DATA,'CellList') && size(DATA.CellList,2) > length(DATA.Expts) * 2 %% old style fr swatches
    cells = unique(DATA.CellList(eid,:,:));
else
    eid = DATA.expid;
    cells = [];
    probelist = DATA.probelist;
end
%eid is now a list of inices to the CellList array
%so expt #s are CellDetails.exptids(eid)
if strfind(DATA.outname,'image.ORBW')
    DATA.plot.flip = ~DATA.plot.flip;
end
DATA.listbycell = 1;
domu = [];

if makemu
    if isfield(DATA,'CellList')
        for p = probelist
            nex = 0;
            nmu = 0;
            mucluster = [];
            for e = 1:length(eid)
                if sum(DATA.CellList(eid(e),p,:)) == 0
                    nex = nex+1;
                    mucluster = 1;
                else
                    mucluster = 0;
                end
                if size(DATA.muCellList,1) >= eid(e)
                    a = find(DATA.muCellList(eid(e),p,:) > 0);
                    if ~isempty(a)
                        mucluster = a;
                        nmu = nmu+1;
                    end
                end
                muclusters(p,e) = mucluster(1);
            end
            if nex >= length(eid)/2 % more htan half expts no cell defined
                domu = [domu p];
            elseif nmu > 0
                domu = [domu p];
            end
            
        end
    else
        domu = probelist;
    end
else
   muclusters = []; 
end

if makemu ==2
    js = length(cells)+1;
else
    js = 1;
end
if oneprobe > 0
    domu = oneprobe;
end
for j = js:(length(cells)+length(domu))
    if j > length(cells)
        mu = domu(j - length(cells));
        DATA.listbycell = 2;
        DATA = cmb.SetProbe(DATA, mu);
        p =mu;
        mucluster = muclusters(p,:);
    else
        DATA = cmb.SetProbe(DATA, cells(j));
        p = cells(j);
    end
    if checkonce
        chk = 0;
    end
    proceed = 1;
    if j > length(cells)
        file = cmb.CombinedName(DATA,DATA.extype,1,'cell',mu,'modifier',modifier);
        file = strrep(file,'.cell','.mu');
    else
        Expt.Header.cellnumber = cells(j);
        file = cmb.CombinedName(DATA,DATA.extype,1,'cell',cells(j),'modifier',modifier);
    end
    if useeachcombineids
        if exist(file)
            load(file);
            combinlist = find(ismember(DATA.expid,Expt.Header.Combineids));;
        else
            fprintf('Cant read %s, so cant determine Expts to use\n',file);
            proceed = 0;
        end
    end
    if proceed
        if p <= size(muclusters,1)
            DATA.muclusters(eid) = muclusters(p,:);
        end
        [Expt, DATA, plotres] = cmb.CombinePlot(DATA, chk,'ids',combinelist);
        Expt.cp = cmb.ExptCP(Expt);
        if savefiles
            ts = now;
            if j > length(cells)
                Expt.Header.cellnumber = 0;
                file = cmb.CombinedName(DATA,DATA.extype,1,'cell',mu,'modifier',modifier);
                file = strrep(file,'.cell','.mu');
            else
                Expt.Header.cellnumber = cells(j);
                file = cmb.CombinedName(DATA,DATA.extype,1,'cell',cells(j),'modifier',modifier);
            end
            c = '';
            nspk = sum([Expt.Trials.count]);
            BackupFile(file,'print');
            save(file,'Expt');
            fprintf('Saved %d spikes (Expts %s) to %s (%s) took %.1f\n',nspk,sprintf(' %d',DATA.combineids),file,DATA.user,mytoc(ts));
        end
        
%Check for any trials that are new in this probe, to make sure AllExpt.Expt has
%all the trials for all probes
        trialids = [AllExpt.Expt.Trials.id];
        [xid, ix] = setdiff([Expt.Trials.id],trialids); %These trials excluded from initial expt.
        if ~isempty(xid)
            fprintf('Adding Ids %s to trial list\n',sprintf(' %d',xid));
            for k = 1:length(xid)
                addedid(end+1) = xid(k);
                AllExpt.Expt.Trials = CopySFields(AllExpt.Expt.Trials,'new',Expt.Trials(ix(k)));
            end
        end
        
        for t = 1:length(Expt.Trials)
            AllExpt.Spikes{j}.trialid(t) = int16(Expt.Trials(t).id);
            AllExpt.Spikes{j}.Trial(t) = int16(Expt.Trials(t).Trial);
            AllExpt.Spikes{j}.Spikes{t} = int16(Expt.Trials(t).Spikes);
            AllExpt.Spikes{j}.OSpikes{t} = int16(Expt.Trials(t).OSpikes);
            AllExpt.Spikes{j}.Ocodes{t} = int8(Expt.Trials(t).Ocodes);
        end
        AllExpt.Header(j).cellnumber = Expt.Header.cellnumber;
        hfields = {'probe' 'probes' 'dropi' 'dips' 'Clusters' 'excludeids' 'Combineids' 'suffixes' 'muforcell'};
        AllExpt.Header = CopySFields(AllExpt.Header,j,Expt.Header,hfields);
        if ismember(Expt.Header.probe,find(DATA.ArrayConfig.badprobes>0))
            AllExpt.Header(j).badprobe = DATA.ArrayConfig.badprobes(Expt.Header.probe);
        else
            AllExpt.Header(j).badprobe = 0;
        end
        
        Expt.plotres = rmfields(plotres,'Data');
        if isfield(DATA,'AllClusters') && ~isfield(Expt.Header,'probesep')
            Expt.Header.probesep = 50;
            AllExpt.Expt.Header.probesep = 50;
        end
        %    Expt.Header.SpkStats = cmb.GetSpkStats(DATA);
        drawnow;
        nspk = sum([Expt.Trials.count]);
        if DATA.listbycell == 2
            file = regexprep(outname,'\.p[0-9]*c1\.',['.mu' num2str(p) 'c1.']);
        elseif DATA.listbycell
            file = regexprep(outname,'\.p[0-9]*c1\.',['.cell' num2str(cells(j)) 'c1.']);
        else
            file = regexprep(outname,'\.p[0-9]*c1\.',['.p' num2str(DATA.probelist(p)) 'c1.']);
        end
        %    file = cmb.CombinedName(DATA,eid,1);
        if DATA.state.nospikes == 0
            BackupFile(file,'print');
            save(file,'Expt','ExptHeader');
            fprintf('Saved %d spikes to %s\n',nspk,file);
        end
        Expts{j} = Expt;
        if isfield(Expt,'suffixes')
            estr = sprintf('%s Suffs %s',sprintf(' %d',Expt.Header.Combineids),sprintf(' %d',Expt.Header.suffixes));
        else
            estr = sprintf(' %d',Expt.Header.Combineids);
        end
        if DATA.logfid > 0
            fprintf(DATA.logfid, '%s,Saved %d spikes (Expts%s) to %s (%s) AllCell%d\n',datestr(now),nspk,estr,file,DATA.user,j);
        end
    end
end
if ~isempty(addedid)
    [a,b] = sort([AllExpt.Expt.Trials.Trial]);
    AllExpt.Expt.Trials = AllExpt.Expt.Trials(b); 
end
ExptHeader = AllExpt.Expt.Header;


if isfield(DATA,'errs')
    AllExpt.errs = DATA.errs;
end
setappdata(DATA.toplevel,'AllExpt',AllExpt)
DATA.Expt = Expt;
DATA.plot.autoclustermode = oldauto;
if DATA.state.online == 0 && dolfp
    cmb.combine('savelfp',DATA);
end
if makemu == 2
    Expts = Expts(js:end);
end
DATA.AllExpts = Expts;
DATA.plot = oldplot;
DATA.listbycell = oldlistbycell;
DATA.state.showspikes = oldspikes;

if guicall
    set(DATA.toplevel,'UserData',DATA);
end
ename = Expt2Name(DATA.Expts{DATA.currentexpt(1)});
%make sure allexptname is something different (i.e. regexprep worked)
if strcmp(DATA.outname, allexptname)
    allexpname = [fileparts(DATA.outname) '/AllCellsTemp.mat'];
    fprintf('Changing name to %s\n',allexptname);
end

if ~strcmp(DATA.outname, allexptname)
    BackupFile(allexptname,'print');
    fprintf('Saving AllExpt to %s\n',allexptname);
    save(allexptname,'AllExpt','ExptHeader')
end
if DATA.state.interactive
    cmb.SetCellChooser(DATA);
    cmb.PlotAllCells(DATA,[],'rates');
end
mytoc(ts);
