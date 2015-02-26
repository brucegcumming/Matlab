function out = combine(name, varargin)
%combine(file)
%cuts clusters and combines expts for Spike2 generated matlab data files
%
%


%
% Check for sameness of size, dx,sf, tf etc (CheckTrl equiv);
% record noclusters
% add checkbox for plotflip
%.
% need to allow other clustering  parameters, esp for pre-dip cels, like
% ruf1907

TOPTAG = 'Combiner'; %default
nextcell = 0;
if length(varargin) & isstruct(varargin{1})
   DATA = varargin{1};
   TOPTAG = DATA.tag.top;
else
    for j = 1:length(varargin)
        if ischar(varargin{j}) && strcmp('Tag',varargin{j})
            TOPTAG = varargin{j+1};
        end
    end
end
it = findobj('Tag',TOPTAG);
if isempty(it)
DATA.expstrs = {'dx' 'xo' 'dp' 'dO' 'dxXce' 'sf' 'dxXId' 'or' 'Op' 'Pp' 'sz' 'jv'...
    'orXob' 'orXme' 'szXob' 'DcXor' 'OpXme' 'PpXme' 'co' 'CtfXip' 'tfXip' 'IB' 'backMov' 'backMovXannTyp' 'IBXannTyp' 'OD'};
DATA.expnames = {'DT', 'XO', 'DP' 'ODX' 'AC' 'SF' 'ABD' 'OT' 'OP' 'PP' 'SZ' 'VE' ...
    'ORBW' 'OXM' 'SZOB' 'DCOR' 'OPM' 'PPM' 'CO' 'CTFIP' 'TFIP' 'IB' 'Movie' 'MovieANN' 'IB' 'OD'};
DATA.expstrs = {DATA.expstrs{:} 'sfXob' 'sfXme' 'szXme' 'dxXdx' 'dOXdO'};
DATA.spkvarnames = GetSpikeVals(DATA,NaN,NaN,NaN);
DATA.expnames = {DATA.expnames{:} 'SFOB' 'SFM' 'SZM' 'DT' 'ODX'};
    DATA.suffs = 'abcdeghijklmnopqrstuvwxyz';
    DATA.state.recut = 0;
    DATA.state.recount = 0;
    DATA.state.plotpsych = 0;
    DATA.state.plotcombined = 0;
    DATA.state.autoplot = 0;
    DATA.state.autolist = 0;
    DATA.state.plotseq  = 0;
    DATA.state.uselfp = 0;
    DATA.tag.top = TOPTAG;
    DATA.tag.options = 'CombineOptions';
    DATA.tag.spikev = 'SpikeV';
    DATA.tag.clusterxy = 'ClusterPlot';
    DATA.tag.emplot = 'Combined EM';
    colors = mycolors;
    DATA.spkcolor{1} = [0.5 0.5 0.5];
    DATA.spkcolor(2:20) = colors(1:19);
    DATA.tag.dataplot = 'CombinerData';
    DATA.datafilename = name;
    DATA.wsc = 1;
    DATA.defaults.fz = 96;
    DATA.state.fixrange = 0;
    DATA.spooling = 0;
    DATA.densityplot = 0;
    DATA.plot.clusterXrange = [0 10];
    DATA.plot.clusterYrange = [0 1.5];
    DATA.plot.clusterZrange = [0 10];
    DATA.cluster{1}.Arange = [6:8];
    DATA.cluster{1}.Brange = [11:20];
    DATA.cluster{1}.Erange = [1:100];
    DATA.clusterArange = [6:8];
    DATA.clusterBrange = [11:20];
    DATA.clusterErange = [1:100];
    DATA.cluster{1}.params = [1 2];
    DATA.plot.autoscale = 1;
    DATA.state.autofit = 0;
    DATA.plot.nodc = 1;
    DATA.test.fastplot = 0;
    DATA.show.sz = 0;
    DATA.show.or = 0;
    DATA.show.jv = 0;
    DATA.show.sf = 0;
    DATA.show.me = 0;
    DATA.show.tf = 0;
    DATA.show.dw = 0;
    DATA.xyfig = 0;
    DATA.optionfig = 0;
    DATA.svfig = 0;
    DATA.plot.showem = 0;
    DATA.plot.showcp = 0;
    DATA.spikelist = 1;
    DATA.plot.dvdt = 0;
    DATA.minplottime = 0.1;
    DATA.plot.clusterX = 1;
    DATA.plot.clusterY = 2;
    DATA.plot.nmin = 0;
    DATA.plot.sdfw = 0;
    DATA.plot.showISI = 0;
    DATA.plot.acov = 0;
    DATA.plot.collapse = 0;
    DATA.plot.flip = 0;
    DATA.plot.bsdelaymax = 5000;
    DATA.plot.setptsize = 0;
    DATA.probenames = {'1' '2'};
    DATA.probelist = [5 4];
    DATA.probe = 5;
    DATA.linesread = 0;
    DATA.currentcluster = 1;
    DATA.state.online = 0;
    DATA.outname = 'tmp.mat';
    DATA.badnames = {};
    DATA.fitvals = [];
    DATA = BuildGUI(DATA);

    firstcall = 1;
else
    firstcall = 0;
    if length(varargin) & isstruct(varargin{1})
        DATA = varargin{1};
        TOPTAG = DATA.tag.top;
    else
        DATA = get(it,'UserData');
    end
    set(DATA.toplevel,'Name','Busy......');
    drawnow;
end

if strncmpi(name,'nextcell',5)
    dir = splitpath(DATA.datafilename,'dir');
    if dir > 99
    newdir = num2str(str2num(dir)+1);
    else
        newdir = num2str(str2num(dir)+1,'%.3d');
    end
    name = strrep(DATA.datafilename,dir,newdir);
end

reindex = 0;
j = 1;
while j <= nargin -1
    if strncmpi(varargin{j},'relist',4)
        reindex = 1;
    end
    j = j+1;
end



if strncmpi(name,'listexps',6)
    eid = get(DATA.clst,'value');
    if eid(1) > 1
        suff = '';
        it = strmatch(DATA.exptypelist{eid(1)},DATA.expstrs,'exact');
        stimname = strrep(DATA.explist{eid(1)},['.' DATA.exptypelist{eid(1)}],'');
        if strfind(stimname,'RC')
            stimname = strrep(stimname,'RC','');
            suff = 'RC';
        end
        if ~isempty(it)
            expname = DATA.expnames{it};
            if strmatch(expname,'CO','exact') & strmatch(stimname,'square')
                expname = 'FLSH';
            end
            if regexp(DATA.datafilename,'Expt[0-9]*.mat') %online file
                DATA.outname = [DATA.datafilename '/c' num2str(DATA.spikelist(1)) '.' stimname '.' DATA.expnames{it} suff '.mat'];
            else
                DATA.outname = CombinedName(DATA,eid(1),DATA.spikelist(1));
            end
            if exist(DATA.outname,'file')
                GetFigure(DATA.tag.dataplot);
                load(DATA.outname);
                Expt.Header.Name = BuildName(Expt.Header.Name);
                DATA.Expt = Expt;
                args = PlotArgs(DATA, Expt);
                DATA.plotres = PlotExpt(Expt,args{:});
                csuffs = [];
                for j = 1:length(Expt.Header.Combined)
                    csuffs = [csuffs num2str(Expt.Header.Combined(j)) ' '];
                end
               title(sprintf('%s %s ed %.2f',splitpath(DATA.outname),csuffs,GetEval(Expt,'ed')));
            end
        else
            DATA.outname = 'tmp.mat';
        end
    else
            DATA.outname = 'tmp.mat';
    end
    DATA = ListSubExpts(DATA,eid);
    set(DATA.saveitem,'string',DATA.outname);
    set(DATA.toplevel,'UserData',DATA);
    if nargout
        out = DATA;
    end
elseif strncmpi(name,'getexpt',6)
    out = DATA.Expt;
    NotBusy(DATA);
    return;
elseif strncmpi(name,'getstate',6)
    out = DATA;
    NotBusy(DATA);
    return;
elseif strncmpi(name,'checklists',6)
    CheckLists(DATA);
elseif strncmpi(name,'combine',6)
    DATA.Expt = CombinePlot(DATA);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'options',5)
    DATA.optionfig = setoptions(DATA,DATA.tag.options);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'PlotISI',5)
    GetFigure('ISI');
    [isis, t, s] = CalcISI(DATA.Expts{DATA.currentexpt}.Trials);
    id = find(isis < 1000)
    hist(isis(id),100);
    [a,b] = sort(isis);
%    DATA.isit = t(b);
    DATA.isis = s(b);
    DATA.ISIpair = 1;
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'showvals',5)
    DATA.showid = setshow(DATA,DATA.tag.options);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'spikes',5)
    PlotSpikes(DATA,varargin{1});
elseif strncmpi(name,'setexp',6)
    ts  = now;
    id = get(DATA.elst,'value');
    GetFigure(DATA.tag.dataplot,'noforce');
    if DATA.state.plotcombined
        CombinePlot(DATA);
        NotBusy(DATA);
        return;
    end

    hs = 'nohold';
    colors = mycolors;
    if strncmpi(name,'setexpplot',8)
        playspk = 0;
    else
        playspk = get(findobj(DATA.toplevel,'Tag','ShowSpikes'),'value');
    end
    DATA.currentexpt = DATA.expid(id(1));
    if DATA.state.online 
        ton = strcmp(get(DATA.timerobj,'Running'),'on');
        if DATA.expid(id(end)) == length(DATA.Expts) || DATA.state.autolist
            if ~ton
                start(DATA.timerobj);
                DATA = get(DATA.toplevel,'UserData');
            end
        elseif ton
            stop(DATA.timerobj);
        end
    end
    if DATA.state.online < 2 & DATA.Expts{DATA.expid(id(1))}.gui.classified == 0 || ...
             (DATA.Expts{DATA.expid(id(1))}.gui.clustertype == 2 && DATA.state.autoplot)
        playspk = 1;
    end
    if playspk
        if DATA.spooling %was looking at a different classification. Re-classif
            for j = id(1:end)
            end
        end
        DATA.spklist = [];
        DATA.spooling = 0;

        for j = id(1:end)
 % If this expt last shown with a temporary cut, revert back to the saved
 % one, unless this is a "spool" showing the temporary cut.
 % or if the classifiction has  not been done at all yet.
            if ~DATA.spooling && ismember(DATA.Expts{DATA.expid(j)}.gui.classified,[0 2])
                DATA = SetExptSpikes(DATA,DATA.expid(id(1)),0);
            else
                [DATA.plot.clusterX , DATA.plot.clusterY, DATA] = GetClusterSpace(DATA, DATA.Expts{DATA.expid(id(1))});        
            end
            DATA = PlaySpikes(DATA, DATA.expid(j));
        end
        set(DATA.toplevel,'UserData',DATA);
        if ~DATA.state.autoplot
            NotBusy(DATA);
            return; % This stops plotting of results. Why? ? only return if not auto
        end
    end
    DATA.spklist = []; %once clusters are cut, this should be in them. Don't let first expt set for all
    GetFigure(DATA.tag.dataplot,'noforce');
    ClearPlot;
    h = [];
    for j = id(1:end)
        if DATA.state.recount | DATA.Expts{DATA.expid(j)}.gui.counted == 0;
            DATA = CountSpikes(DATA,DATA.expid(j));
        end
        te = now;
            args = PlotArgs(DATA,DATA.Expts{DATA.expid(j)});
        res = PlotExpt(DATA.Expts{DATA.expid(j)},hs,'forcecolor',colors{j},args{:},'legendpos',7);
%        fprintf('%.2f ',(now-te)*60*60*24);
        sp  = strfind(DATA.subexplist{j},' ');
        if sp
            names{j} = DATA.subexplist{j}(1:sp(1)-1);
        else
            names{j} = DATA.subexplist{j};
        end
        fn = fields(DATA.show);
        for k = 1:length(fn)
            if DATA.show.(fn{k})
                names{j} = [names{j} sprintf(' %s=%.2f',fn{k},GetEval(DATA.Expts{DATA.expid(j)},fn{k}))];
            end
        end
        if  ~isempty(res) && isfield(res,'handles') && ishandle(res(1).handles(1))
            h(j) = res(1).handles(1);
        end
        hs = 'Hold';
    end
    if ~isempty(h)
        mylegend(h,names);
    end
    t = get(get(gca,'title'),'String');
    title([t 'Cl' sprintf(' %d',WhichClusters(DATA.toplevel))]);
 %   fprintf('Took %.2f\n',(now-ts)*60*60*24);
    set(DATA.toplevel,'UserData',DATA);
    
elseif strncmpi(name,'Close',4)
    if isfield(DATA,'timerobj')
        stop(DATA.timerobj);
    end
    if isfield(DATA,'tag')
        names = fieldnames(DATA.tag);
        for j = 1:length(names)
            CloseTag(DATA.tag.(names{j}));
        end
    end
    return;
elseif strncmpi(name,'save',4)
    DATA.outname = get(DATA.saveitem,'string');
    Expt = DATA.Expt;
    Expt.Header.CombineDate = now;
    if isfield(DATA,'cluster')
        Expt.Header.Cluster = DATA.cluster;
    end
    Expt.Header.clist = DATA.spikelist;
    save(DATA.outname,'Expt');
    cfile = CombinerLst(DATA);
    if exist(cfile,'file') == 2
        load(cfile);
    else
        combines.name = {};
    end
    id = strmatch(splitpath(DATA.outname),combines.name,'exact');
    if isempty(id)
        id = length(combines.name)+1;
    end
    combines.name{id} = splitpath(DATA.outname);
    combines.lastdate(id) = now;
    combines.list{id} = Expt.Header.Combined;
    combines.starts{id} = Expt.Header.BlockStart;
    combines.trials(id,:) = [Expt.Trials(1).Trial Expt.Trials(end).Trial];
    DATA.combines = combines;
    %tic;
    DATA = ListExpts(DATA,DATA.Expts);
   % toc
    set(DATA.toplevel,'UserData',DATA);
  %  combine('listexpts');
  %  save(cfile,'combines');
elseif strncmpi(name,'store',5)
    set(DATA.toplevel,'UserData',varargin{1});
elseif strncmpi(name,'relist',6)
    DATA = ReadDir(DATA, DATA.name);
    set(DATA.toplevel,'UserData',DATA);
elseif exist(name,'dir') % do dir before file, since dirs pass exist(name,'file')
    DATA.Expts = {};
    args = {};
    if reindex
        args = {args{:}, 'relist'};
    end
    DATA.state.onine = 1;
    DATA = ReadDir(DATA, name, args{:});
    DATA = LoadClusters(DATA,ClusterFile(DATA));
    DATA.lastread = 0;
    DATA.lastsize = 0;
    if ~isfield(DATA,'timerobj')
            DATA.timerobj = timer('TimerFcn',@timerfna, 'Period', 1.0, 'ExecutionMode','FixedSpacing');
            DATA.lastread = 0;
    end
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'newfile',6)
    DATA.datafilename = get(findobj(DATA.toplevel, 'Tag','FileName'),'string');
    combine(DATA.datafilename);
elseif strncmpi(name,'rfupdate',6)
    if isfield(DATA.fitvals,'Pp') & isfield(DATA.fitvals,'Op')
        if abs(DATA.fitvals.OpRo-DATA.fitvals.PpRo) > 2
            message('Ro mismatch');
        end
        xy = op2xy([DATA.fitvals.Op DATA.fitvals.Pp],DATA.fitvals.OpRo);
        uflname = strrep(DATA.datafilename,'.mat','.ufl');
        fid = fopen(uflname,'a');
        if fid > 0
        fprintf(fid,'FitRF %.2f %.2f %.2f %.2f\n',xy(1),xy(2),DATA.fitvals.Opw,DATA.fitvals.Ppw);
        fclose(fid);
        end
    end
elseif strncmpi(name,'winfront',6)
    f = fields(DATA.tag);
    for j = 1:length(f)
        it = findobj('Tag',DATA.tag.(f{j}));
        if ~isempty(it)
            figure(it);
        end
    end
elseif exist(name,'file')
    args = {};
    if reindex
        args = {args{:}, 'relist'};
    end
    DATA.datafilename = name;
    DATA.state.online = 0;
    set(findobj(DATA.toplevel, 'Tag','FileName'),'string',name);
    if strfind(name,'.txt') %% and online text file
        Expts = ReadOnlineTxt(name, args{:});
        Expts = CountTxtSpikes(Expts,DATA.probe,DATA.spikelist);
        DATA.state.online = 2;
        DATA.lastread = 0;
        DATA.lastsize = 0;
        if ~isfield(DATA,'timerobj')
            DATA.timerobj = timer('TimerFcn',@timerfn, 'Period', 1.0, 'ExecutionMode','FixedSpacing');
            DATA.lastread = 0;
        end
    else
        [Trials, Expts, All] = APlaySpkFile(name, args{:});
    end
    DATA = ListExpts(DATA,Expts);
    Trialids = [];
    for nexp = 1:length(Expts)
        newt = [Expts{nexp}.Trials.Trial];
        Trialids = [Trialids newt];
        Expts{nexp}.gui.classified = 0;
        Expts{nexp}.gui.counted = 0;
        Expts{nexp}.gui.clustertype = 0;
        if isfield(Expts{nexp}.Header,'Nam1e')
            Expts{nexp}.Header.Name = Expts{nexp}.Header.Nam1e;
        end
    end
    if isempty(Expts)
        NotBusy(DATA);
        return;
    end
    DATA.Expts = Expts;
    if DATA.state.online
        DATA.AllData = [];
    else
    DATA.AllData = All;
    DATA.AllData.SpikeIdx = Trials.Spkid;
    DATA.AllData.Trialids = Trialids;

    
% load saved cluster params
    
    DATA = LoadOnlineClusters(DATA, ClusterFile(DATA,'getonline'));
    DATA = LoadClusters(DATA, ClusterFile(DATA));
    cfile = CombinerLst(DATA);
    if exist(cfile,'file')
        load(cfile);
        DATA.combines = combines;
    end
    end

    set(DATA.toplevel,'UserData',DATA);
    if strfind(name,'ic-169')
        CheckSaccades(DATA.Expts,'Z:/smr/icarus/169/ic169.em.mat');
    end
elseif regexp(name,'ruf[0-9][0-9][0-9]') %if its not a file name, try building path
    dnum = strrep(name,'ruf','');
    set(DATA.toplevel,'UserData',DATA);
    name = ['/data/rufus/' dnum '/' name '.mat'];
    if exist(name,'file')
        combine(name,'Tag',DATA.tag.top);
    end
elseif firstcall
    fprintf('No file or directory %s\n',name);
    close(DATA.toplevel);
    return;
end
SetGui(DATA);
NotBusy(DATA);

function NotBusy(DATA)

[a,b,c] = fileparts(DATA.datafilename);
set(DATA.toplevel,'Name',[DATA.tag.top ':' b]);

function Expts = CountTxtSpikes(Expts, probe, cl)

if length(cl) > 2
    cl = cl(1:2);
end
for j = 1:length(Expts)
    if isfield(Expts{j},'Trials')
    for k = 1:length(Expts{j}.Trials)

        if size(Expts{j}.Trials(k).AllSpikes,2) < cl(1)+1
            Expts{j}.Trials(k).Spikes = [];
        elseif length(cl) > 1 && size(Expts{j}.Trials(k).AllSpikes,2) > cl(end)
            Expts{j}.Trials(k).Spikes = union(Expts{j}.Trials(k).AllSpikes{probe,cl+1});
        else
            Expts{j}.Trials(k).Spikes = Expts{j}.Trials(k).AllSpikes{probe,cl(1)+1};
        end
    end
    else
        fprintf('Expt %d no trials\n');
    end
    Expts{j}.gui.counted = 1;
end
        
function DATA = LoadClusters(DATA, cfile)
    if exist(cfile,'file')
        excludelist = [];
        clid = [];
        classified = 0;
        load(cfile);
        DATA.state.recut = 1;
        firstspk = length(DATA.AllData.Spikes.times);
        lastspk = 0;
        if size(clid,1) == size(DATA.AllData.Spikes.codes,1);
            DATA.AllData.Spikes.codes(:,2) = clid;
            last = max(find(clid > 0));
            maxclasst = DATA.AllData.Spikes.times(last);
            classified = 1;
        end
        for j = 1:min([length(excludelist) length(DATA.Expts)]);
            lst = -abs([DATA.Expts{j}.Trials(excludelist{j}).Trial]);
            if lst
                for k = 1:length(lst)
                    [DATA.Expts{j}.Trials(excludelist{j}(k)).Trial] = lst(k);
                end
            end
        end
        for j = 1:min([length(AllClusters) length(DATA.Expts)]);
            if ~isempty(AllClusters{j})
                DATA.Expts{j}.Cluster = AllClusters{j}.Cluster;
                if exist('clustertypes','var') %%saved type
                    DATA.Expts{j}.gui.clustertype = clustertypes(j);
                else
                    DATA.Expts{j}.gui.clustertype = 1;
                end
        firstspk = length(DATA.AllData.Spikes.times);
        lastspk = 0;
        nset = 0;
                for k = 1:length(DATA.Expts{j}.Cluster)
                    if isempty(DATA.Expts{j}.Cluster{k}) | ~isfield(DATA.Expts{j}.Cluster{k},'firstspk')
                        DATA.Expts{j}.Cluster{k}.firstspk = NaN;
                        DATA.Expts{j}.Cluster{k}.lastspk = NaN;
                    else
                        nset = nset + 1;
                    end
                    if DATA.Expts{j}.Cluster{k}.firstspk < firstspk
                        firstspk = DATA.Expts{j}.Cluster{k}.firstspk;
                    end
                    if DATA.Expts{j}.Cluster{k}.lastspk > lastspk
                        lastspk = DATA.Expts{j}.Cluster{k}.lastspk;
                    end
                    if DATA.Expts{j}.gui.clustertype == 2 % online cluster - can't use spk coun
                           DATA.Expts{j}.Cluster{k}.firstspk = NaN;
                           DATA.Expts{j}.Cluster{k}.lastspk = NaN;
                    end
                    if ~isfield(DATA.Expts{j}.Cluster{k},'params')
                        DATA.Expts{j}.Cluster{k}.params = [1 2];
                    end                    
                    if ~isfield(DATA.Expts{j}.Cluster{k},'Arange')
                        DATA.Expts{j}.Cluster{k}.Arange = DATA.cluster{1}.Arange;
                        DATA.Expts{j}.Cluster{k}.Brange = DATA.cluster{1}.Brange;
                        DATA.Expts{j}.Cluster{k}.Erange = DATA.cluster{1}.Erange;
                    end                    
                end
                if lastspk < firstspk %nevet found these in cluster file
                    firstspk = NaN;
                    lastspk = NaN;
                end
                if nset == 0 % actively set no clusters
                    DATA.Expts{j}.gui.classified = 1;
                    DATA.Expts{j}.gui.ncluster = 0;
                elseif classified && lastspk > firstspk && max(DATA.AllData.Spikes.codes(firstspk:lastspk,2)) > 0 
                    DATA.Expts{j}.gui.classified = classified;
                    DATA.Expts{j}.gui.ncluster = nset;
                else
                    DATA.Expts{j}.gui.classified = 0;
                    DATA.Expts{j}.gui.ncluster = 0;
                end
            end
        end
    end

function DATA = LoadOnlineClusters(DATA, cfile)
    if exist(cfile,'file')
        load(cfile);
        DATA.state.recut = 1;
        for j = 1:length(AllClusters);
            if ~isempty(AllClusters{j})
                for k = 1:length(DATA.Expts)
                    if isfield(AllClusters{j},'ids') & ...
                            AllClusters{j}.ids(1) > min([DATA.Expts{k}.Trials.id]) -2 & ...
                            AllClusters{j}.ids(2) < max([DATA.Expts{k}.Trials.id]) +2
                        
                        DATA.Expts{k}.Cluster = AllClusters{j}.Cluster;
                        DATA.Expts{k}.gui.clustertype = 2;
%
%Fill in fields that might be missing from old cluster files
%remove firstpk, lastspk refs from the online file - these spike id #s will
%not necessarily match the saved fil
                        for m = 1:length(DATA.Expts{k}.Cluster)
                            DATA.Expts{k}.Cluster{m}.firstspk = NaN;
                            DATA.Expts{k}.Cluster{m}.lastspk = NaN;
                            if ~isfield(DATA.Expts{k}.Cluster{m},'params')
                                DATA.Expts{k}.Cluster{m}.params = [1 2];
                            end
                        end
                        DATA.Expts{k}.gui.classified = 0; %loaded def, but not classified
                    end
                end
            end
        end
    end


function DATA = CheckCombine(DATA)
    id = get(DATA.elst,'value');
    for j = 1:length(id)
        stims(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'st');
        bws(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'ob');
        mes(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'me');
        tfs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'tf');
        sfs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'sf');
        ets{j} = DATA.Expts{DATA.expid(id(j))}.Stimvals.et;
        e2s{j} = DATA.Expts{DATA.expid(id(j))}.Stimvals.e2;
    end
    for j = 1:length(e2s)
            if isempty(e2s{j}) 
            e2s{j} = 'e0';
            end
            if isempty(ets{j}) 
            ets{j} = 'e0';
            end
    end
    blank = strmatch('none',unique(stims));
    nstims = length(unique(stims)) - length(blank);
    blank = strmatch('e0',unique(e2s));
    e2lst = unique(e2s);
    e1lst = unique(ets);
    ne2 = length(e2lst) - length(blank);
    if bws > 1 & ne2 == 0 
        for j = id(1:end)
            DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},'ob');
            DATA.Expts{DATA.expid(j)}.Stimvals.e2 = 'ob';
        end
    elseif ne2 == 1 & strmatch('me',e2lst)
        for j = id(1:end)
            DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},'me');
            DATA.Expts{DATA.expid(j)}.Stimvals.e2 = 'me';
        end        
    end
    if length(unique(sfs)) >1 & isempty(strmatch('sf',{e1lst{:} e2lst{:}}))
        questdlg([sprintf('%d SFS in %s',length(unique(sfs)),DATA.outname) sprintf('%.2f ',unique(sfs))],'test','OK','OK');
    end
    

function Expt = CombinePlot(DATA)
    
    if DATA.state.online < 2
        DATA = CheckCombine(DATA);
    end
    id = get(DATA.elst,'value');
    tfields = {};
    for j = id(1:end)
        tfields = union(tfields,fields(DATA.Expts{DATA.expid(j)}.Trials));
    end
    j = DATA.expid(id(end));
    Expt.Header = DATA.Expts{j}.Header;
    Trials = [];
    Pulses = [];
    nb = 1;
    Clusters = {};
    for j = id(1:end)
        BlockStart(nb) = DATA.Expts{DATA.expid(j)}.Trials(1).Trial;
        Header = DATA.Expts{DATA.expid(j)}.Header;
        DATA = CountSpikes(DATA,DATA.expid(j)); %% recount in case clustere # changed
        if isfield(DATA.Expts{DATA.expid(j)},'Cluster')
            Clusters{nb} = DATA.Expts{DATA.expid(j)}.Cluster;
        end
%tfields is the combined fields of all expts (above). So any differences
%mean it is missing from the current Expt;
        newf = setdiff(tfields, fields(DATA.Expts{DATA.expid(j)}.Trials));
        newtrials = DATA.Expts{DATA.expid(j)}.Trials;
        for k = 1:length(newf)
            if strcmp(newf{k},'FalseStart')
            [newtrials.(newf{k})] = deal(0);
            else
            [newtrials.(newf{k})] = deal([]);
            end
        end
        Trials = [Trials newtrials];
        if isfield(DATA.Expts{DATA.expid(j)},'ExcludeId')
        end
        if isfield(DATA.Expts{DATA.expid(j)},'Pulses')
            Pulses = [Pulses DATA.Expts{DATA.expid(j)}.Pulses];
        end
        if Header.trange(1) < Expt.Header.trange(1)
            Expt.Header.trange(1) = Header.trange(1);
        end
        if Header.trange(2) > Expt.Header.trange(2)
            Expt.Header.trange(2) = Header.trange(2);
        end
        nb = nb+1;
    end
    j = DATA.expid(id(end));
    Expt.Stimvals = DATA.Expts{j}.Stimvals;
    Expt.Trials = Trials;
    Expt.Header.Name = BuildName(DATA.Expts{j}.Header.Name);
    Expt.Header.BlockStart = BlockStart;
    Expt.Header.Clusters = Clusters;
    Expt.Header.Combined = id;
    Expt.Header.Spikelist = DATA.spikelist;
    if isfield(Expt.Trials,'FalseStart')
        id = find([Expt.Trials.FalseStart] > 0);
        if ~isempty(id)
            for j = 1:length(id)
                Expt.Trials(id(j)).Start = Expt.Trials(id(j)).Start - Expt.Trials(id(j)).FalseStart;
            end
            msgbox(sprintf('%d Trials with long delays',length(id)));
        end
    end
    if ~isempty(Pulses)
        Expt.Pulses = Pulses;
    end
    if ~isfield(Expt.Header,'Spike2Version')
        Expt.Header.Spike2Version = 1.0;
    end
    
    if DATA.state.uselfp & length(id) > 0  %reload LFP to match lengths etc
        Expt = LoadSpike2LFP(Expt,'reload');
    end
    id = strmatch(DATA.Expts{j}.Header.expname,DATA.expstrs,'exact');
    PlotCombined(DATA, Expt);
  
function [Expt, res] = PlotCombined(DATA, Expt)
    
 if isfield(Expt,'ExcludeCluster')
    for j = 1:length(Expt.ExcludeCluster{1})
        id = find([Expt.Trials.Trial] == Expt.ExcludeCluster{1}(j));
        if id
            Expt.Trials(id).Trial = -Expt.Trials(id).Trial;
        end
    end
 end
    GetFigure(DATA.tag.dataplot);
    args = PlotArgs(DATA,Expt);
    res = PlotExpt(Expt,args{:});
    if DATA.state.plotseq == 4 & isfield(res,'fig')
        set(res(1).fig, 'WindowButtonDownFcn',@FigButtonPressed);
        set(res(1).fig, 'WindowButtonUpFcn',@FigButtonReleased);
        dat = get(res(1).fig,'UserData');
        dat.parentfigtag = DATA.tag.top;
        set(res(1).fig,'UserData',dat);
    end
    t = get(get(gca,'title'),'String');
    title([t 'Cl' sprintf(' %d',WhichClusters(DATA.toplevel))]);

    if DATA.plot.showem
        GetFigure(DATA.tag.emplot);
        Expt = LoadEmData(Expt);
        PlotExptEM(Expt);
    end
    if DATA.state.uselfp
        GetFigure('LFP');
        hold off;        
        CalcLFPPulse(Expt,DATA.AllData,'plot');
    end
        
%    save(DATA.outname,'Expt'); %use save to do the saving....

    
   
function DATA = ReadDir(DATA, name, varargin)  %% Online Plots

d = dir(name);
reindex = 0;
args = {};
j = 1;
while j <= nargin-2
        if strncmpi(varargin{j},'relist',3)
            reindex =1;
            args = {args{:} 'relist'};
        end
    j = j+1;
end
    expnames = {};
    if isempty(DATA.Expts)
        nexp = 1;
        SpkId = [];
        Spikes = [];
        Trialids = [];
    else
        for j=1:length(DATA.Expts)
            expnames{j} = splitpath(DATA.Expts{j}.Header.Name);
 %           DATA.Expts{j}.gui.classified = 0;
 % don't reset this here - resets clusters for files already read. Reset
 % it below only when a file is read.
        end
        nexp = j+1;
        SpkId = DATA.AllData.SpikeIdx;
        Spikes = DATA.AllData.Spikes;
        Trialids = DATA.AllData.Trialids;
        All = DATA.AllData;
        Expts = DATA.Expts;
    end
    DATA.defaults.starttrial = 1;
    lastn = nexp;
    nbad = length(DATA.badnames)+1;
    badidx = 1:nbad-1;
    for j = 1:length(d)
        if regexp(d(j).name,'Expt[0-9]*.mat') & ...
                d(j).bytes > 128 & .....
                isempty(strfind(d(j).name,'idx.mat')) & ...    %exclude the .idx file
                d(j).datenum < now-(1/(24 * 60 * 60)) & ... %at least  1 sec old
                isempty(strmatch(d(j).name,{expnames{:} DATA.badnames{badidx}}))     %dont read if we already have
            [trls, exps, All] = APlaySpkFile([name '/' d(j).name],'Defaults',DATA.defaults,'online',args{:});
            if isempty(exps)
                fprintf('%s No expts\n',d(j).name);
                DATA.badnames{nbad} = d(j).name;
                nbad = nbad+1;
            else
                Expts{nexp} = exps{1};
                Expts{nexp}.gui.classified = 0;
                Expts{nexp}.gui.counted = 0;
                Expts{nexp}.gui.clustertype = 0;
                SpkId = [SpkId; trls.Spkid];
                newt = [Expts{nexp}.Trials.Trial];
                Trialids = [Trialids newt];
                if nexp > 1
                    Spikes.values = [Spikes.values; All.Spikes.values];
                    Spikes.dVdt = [Spikes.dVdt; All.Spikes.dVdt];
                    Spikes.codes = [Spikes.codes; All.Spikes.codes];
                    Spikes.times = [Spikes.times; All.Spikes.times];
                else
                    Spikes = All.Spikes;
                end

                nexp = nexp+1;
                DATA.defaults.starttrial = 1+ exps{1}.Trials(end).Trial;
            end
        end
    end
    if nexp == 1
        questdlg(sprintf('No expts in %s',name),'test','OK','OK');
    else
        DATA = ListExpts(DATA,Expts);
        DATA.Expts = Expts;
        DATA.AllData = All;
        DATA.AllData.Spikes = Spikes;
        DATA.AllData.SpikeIdx = SpkId;
        DATA.AllData.Trialids = Trialids;
        if DATA.state.autoplot && DATA.state.online == 1
            for j = lastn:nexp-1
              DATA = SetExptSpikes(DATA, j, 0);
            end
        end
    end
    DATA.name = name;
    DATA.state.online = 1;
    

    
function PlotSpike(DATA, ispk)

    set(0,'CurrentFigure',DATA.svfig);
    j = DATA.AllData.Spikes.codes(ispk,2)+1;
    if DATA.plot.dvdt == 2
          set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.dVdt(ispk,:),'Xdata',DATA.AllData.Spikes.values(ispk,2:end));
    elseif DATA.plot.dvdt
          set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.dVdt(ispk,:),'Xdata',[1:size(DATA.AllData.Spikes.dVdt,2)]);
    else
        set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.values(ispk,:),'Xdata',[1:size(DATA.AllData.Spikes.values,2)]);
    end
    title(sprintf('%d: Cl %d at %.3f',ispk,j-1,DATA.AllData.Spikes.times(ispk)/10000)); 
    drawnow;

function PlayNextSpike(a,b)
%DATA = combine('getstate');
DATA = GetDataFromFig(a);
if DATA.ISIpair
    PlotISIPair(DATA,DATA.ISIpair);
    DATA.ISIpair = DATA.ISIpair+1;
    set(DATA.toplevel,'UserData',DATA);
    return;
end
set(0,'CurrentFigure',DATA.svfig);
PlotSpike(DATA,DATA.currentspike);
set(0,'CurrentFigure',DATA.xyfig);
PlotSpikeXY(DATA,DATA.currentspike,DATA.spkcolor{DATA.AllData.Spikes.codes(DATA.currentspike,2)+1});
DATA.currentspike = DATA.currentspike+1;
set(DATA.toplevel,'UserData',DATA);

function PlotISIPair(DATA, pair)

GetFigure('SpikeV');
hold off;
spk = find(DATA.AllData.Spikes.times >= floor(DATA.isis(pair)));
spk = spk(1);
tsmp = [1:size(DATA.AllData.Spikes.values,2)] * DATA.AllData.Spikes.interval * 10000;
plot(tsmp,DATA.AllData.Spikes.values(spk-1,:));
[a,ta] = max(DATA.AllData.Spikes.values(spk-1,:));
[b,tb] = max(DATA.AllData.Spikes.values(spk,:));
    hold on;
isi = diff(DATA.AllData.Spikes.times([spk-1 spk]));
times = isi+tsmp+20000*DATA.AllData.Spikes.interval;
isi = isi + (tb-ta) * DATA.AllData.Spikes.interval * 10000;
plot(times,DATA.AllData.Spikes.values(spk,:));
title(sprintf('Time %.0f ISI %.1f',DATA.AllData.Spikes.times(spk-1),isi));

function CutTrial(a,b)
%DATA = combine('getstate');
DATA = GetDataFromFig(a);

t = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial;
DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial = -abs(t);
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'string',sprintf('%d|',[DATA.Expts{DATA.currentexpt}.Trials.Trial]),'value',1);
set(DATA.toplevel,'UserData',DATA);
PlayOneTrial(DATA,a,1);
        
function PlayLastTrial(a, b)
DATA = GetDataFromFig(a);
PlayOneTrial(DATA, a,-1);

function PlayNextTrial(a, b)
DATA = GetDataFromFig(a);
PlayOneTrial(DATA,a,1);

function SelectTrial(a, b)
DATA = GetDataFromFig(a);
c = get(findobj('Tag','ChooseTrial'),'value');
PlayOneTrial(DATA,c,0);

function PlayOneTrial(DATA, a, b)
%DATA = combine('getstate');
expid = DATA.currentexpt;

if b < 0  & DATA.currenttrial > 1 %step back
    DATA.currenttrial = DATA.currenttrial-1;
    while DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial < 0 & DATA.currenttrial > 1 
        DATA.currenttrial = DATA.currenttrial-1;
    end
elseif b > 0 & DATA.currenttrial < length(DATA.Expts{DATA.currentexpt}.Trials) %step back
    DATA.currenttrial = DATA.currenttrial+1;
    while DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial < 0 &  ...
            DATA.currenttrial < length(DATA.Expts{DATA.currentexpt}.Trials) 
        DATA.currenttrial = DATA.currenttrial+1;
    end
elseif b == 0 % step to a trial on list
    DATA.currenttrial = find([DATA.Expts{DATA.currentexpt}.Trials.Trial] == a);    
    DATA.currenttrial = a;    
end
if DATA.currenttrial <= 0
    DATA.currenttrial = 1;
end
Trial = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial);
if Trial.Trial < 0 && b == 0 %%manually select -Trial = include again
    DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial = abs(Trial.Trial);
    it = findobj(DATA.svfig,'Tag','ChooseTrial');
    set(it,'string',sprintf('%d|',[DATA.Expts{expid}.Trials.Trial]),'value',1);
end
itrial = find(DATA.AllData.Trialids == abs(Trial.Trial));
set(0,'CurrentFigure',DATA.svfig);
hold off;
%DATA = PlotTrialSpikes(DATA, itrial, mycolors, DATA.clusters);
if DATA.state.recut
    nc = length(DATA.cluster)+1;
else
    nc = DATA.s2clusters;
end
Trial.ed = GetEVal(DATA.Expts{DATA.currentexpt},'ed',DATA.currenttrial);
DATA = APlotTrialSpikes(DATA, [Trial.Start(1) Trial.End(end)], mycolors, nc, 0,'Trial',Trial);

if DATA.state.uselfp
    PlotLFP(DATA.state,Trial,DATA.Expts{expid}.Header.LFPsamplerate);
end
set(DATA.toplevel,'UserData',DATA);
set(0,'CurrentFigure',DATA.svfig);
tid = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).id;
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'value',DATA.currenttrial);
    
function DATA = PlotTrialSpikes(DATA, itrial, colors, clusters)
    spka = DATA.AllData.SpikeIdx(itrial,1);
    spkb = DATA.AllData.SpikeIdx(itrial,2);
    Spks = DATA.AllData.Spikes;
    if spka
        for spk = spka:spkb;
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            hold on;
            energy  = sum(diff(adc).^2);
            DATA.Spikes.energy(spk)= energy;
            svar(spk) = var(adc);
            DATA.Spikes.vw(spk) = svar(spk)/energy;
        end
        title(sprintf('Trial %d (id%d)',DATA.AllData.Trialids(itrial),...
            DATA.AllData.Trialids(itrial)));

        drawnow;


        set(0,'CurrentFigure',DATA.xyfig);
        for j = 1:length(clusters)
            sp = intersect(clusters{j},[spka:spkb]);
            plot(DATA.Spikes.energy(sp),DATA.Spikes.vw(sp),'.','color',colors{j});
        end
    end

function [x,DATA] = GetSpikeVals(DATA, ispk, type, recalc)

SPKENERGY=1;
SPKVARE=2;
SPKMAXRATE = 3;
SPKPREMINRATE = 4;
SPKPEAKTIME = 5;
SPKPREMIN = 6;
SPKPEAK = 7;
SPKVAR = 8;
SPKSYMMETRY = 9;
SPKCENTROID = 10;
SPKCENTROID = 10;
SPKMINRATE = 11;
SPKMAXRATEA = 12;
SPKMINRATEA = 13;
SPKMEANRATEA = 14;
SPKMAXRATEB = 15;
SPKMINRATEB = 16;
SPKMEANRATEB = 17;
SPKMEANA = 18;
SPKMEANB = 19;
SPKMAXA = 20;
SPKMAXB = 21;
SPKMINA = 22;
SPKMINB = 23;
SPKVARTESTA= 24;
SPKVARTESTB= 25;

if isnan(ispk) %% return names
    x = {'Energy' 'Var/Energy' 'MaxRate' 'PreMinRate' 'PeakTime' 'PreMin' 'Peak' 'Var' 'Symmetry' 'Centroid' 'Minrate' 'Maxrate(A)' ... 
    'Minrate(A)' 'Meanrate(A)' 'Maxrate(B)' 'Minrate(B)' 'Meanrate(B)' 'Mean(A)'...
    'Mean(B)' 'Max(A)' 'Min(A)' 'Max(B)' 'Min(B)' 'test1' 'test2'};
  return;
end



if ~isfield(DATA.cluster{DATA.currentcluster},'Arange') & isfield(DATA.cluster{1},'Arange')
    DATA.cluster{DATA.currentcluster}.Arange = DATA.cluster{1}.Arange;
    DATA.cluster{DATA.currentcluster}.Brange = DATA.cluster{1}.Brange;
    DATA.cluster{DATA.currentcluster}.Erange = DATA.cluster{1}.Erange;
end
arange = DATA.clusterArange;
brange = DATA.clusterBrange;
splen = size(DATA.AllData.Spikes.values,2);
erange = intersect(DATA.clusterErange,1:splen-1);    
    

if type == SPKENERGY
    if recalc
        x  = sum(DATA.AllData.Spikes.dVdt(ispk,erange).^2,2);
        DATA.Spikes.energy(ispk)= x;
    else
        x = DATA.Spikes.energy(ispk);
    end
elseif type == SPKVARE
    if recalc
        x = var(DATA.AllData.Spikes.values(ispk,erange)');
        x = x./DATA.Spikes.energy(ispk);
        DATA.Spikes.vw(ispk) = x;
    else
        x = DATA.Spikes.vw(ispk);
    end
elseif type == SPKMAXRATE
    x = max(DATA.AllData.Spikes.dVdt(ispk,:)');
elseif type == SPKPREMINRATE
    [x,t] = max(DATA.AllData.Spikes.dVdt(ispk,:)');
    for j = 1:length(t)
        x(j) = min(DATA.AllData.Spikes.dVdt(ispk(j),1:t(j)));
    end
elseif type == SPKMAXRATEA
    x = max(DATA.AllData.Spikes.dVdt(ispk,arange)');
elseif type == SPKMAXRATEB
    x = max(DATA.AllData.Spikes.dVdt(ispk,brange)');
elseif type == SPKMINRATE
    x = min(DATA.AllData.Spikes.dVdt(ispk,:)');
elseif type == SPKMINRATEA
    x = min(DATA.AllData.Spikes.dVdt(ispk,arange)');
elseif type == SPKMINRATEB
    x = min(DATA.AllData.Spikes.dVdt(ispk,brange)');
elseif type == SPKMEANA
    x = mean(DATA.AllData.Spikes.values(ispk,arange)');
elseif type == SPKMEANB
    x = mean(DATA.AllData.Spikes.values(ispk,brange)');
elseif type == SPKMEANRATEA
    x = mean(DATA.AllData.Spikes.dVdt(ispk,arange)');
elseif type == SPKMEANRATEB
    x = mean(DATA.AllData.Spikes.dVdt(ispk,brange)');
elseif type == SPKMAXA
    x = max(DATA.AllData.Spikes.values(ispk,arange)');
elseif type == SPKMAXB
    x = max(DATA.AllData.Spikes.values(ispk,brange)');
elseif type == SPKMINA
    x = min(DATA.AllData.Spikes.values(ispk,arange)');
elseif type == SPKMINB
    x = min(DATA.AllData.Spikes.values(ispk,brange)');
elseif type == SPKPREMIN
    [x,t] = max(DATA.AllData.Spikes.dVdt(ispk,:)');
    for j = 1:length(t)
        x(j) = min(DATA.AllData.Spikes.values(ispk(j),1:t(j)));
    end
elseif type == SPKPEAKTIME
    [x,t] = max(DATA.AllData.Spikes.dVdt(ispk,:)');
    zc = diff(sign(DATA.AllData.Spikes.dVdt(ispk,:)')); 
    len = size(zc,1);
    for j = 1:length(t)
        tm = find(zc(:,j) < 0 & [1:len]' >= t(j));
        if isempty(tm)
            x(j) = 0;
        else
            x(j) = tm(1); % first zero crossing after peakrate
        end
    end
elseif type == SPKPEAK
    [x,t] = max(DATA.AllData.Spikes.values(ispk,:)');
elseif type == SPKVAR
    x = var(DATA.AllData.Spikes.values(ispk,:)');
elseif type == SPKSYMMETRY
  len =  size(DATA.AllData.Spikes.dVdt,2);
  cn  = round(sum((DATA.AllData.Spikes.dVdt(ispk,:).^2)*[1:len]',2)./sum(DATA.AllData.Spikes.dVdt(ispk,:).^2,2));
  for j = 1:length(ispk)
  z(j) = mean(DATA.AllData.Spikes.values(ispk(j),max([cn(j)-5 1]):cn(j)),2);
  y(j) = mean(DATA.AllData.Spikes.values(ispk(j),cn(j):min([cn(j)+5 len])),2);
  end
  x = atan2(y,z);
elseif type == SPKCENTROID
  len =  size(DATA.AllData.Spikes.dVdt,2);
  x  = sum((DATA.AllData.Spikes.dVdt(ispk,:).^2)*[1:len]',2)./sum(DATA.AllData.Spikes.dVdt(ispk,:).^2,2);
end

function DATA = PlotSpikes(DATA, ispk)

    classify = 1;
    ctype = 2;
    [cx, DATA] = GetSpikeVals(DATA,ispk, DATA.plot.clusterX, classify);
    DATA.Spikes.cx(ispk) = cx;
    %      [cy, DATA] = GetSpikeVals(DATA,ispk, SPKVARE, classify);
    [cy, DATA] = GetSpikeVals(DATA,ispk, DATA.plot.clusterY, classify);
    DATA.Spikes.cy(ispk) = cy;
    DATA.currentspike = ispk(1);
    DATA = SetSpkCodes(DATA,ispk,0);
    nc =  length(DATA.cluster)+1;
    DrawSpikeWaves(DATA, ispk, nc, 2);
    
    set(0,'CurrentFigure',DATA.xyfig);
    for j = 0:nc
        sp = find(DATA.AllData.Spikes.codes(ispk, ctype) == j);
        plot(cx(sp),cy(sp),...
            '.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
        hold on; %% need this when called from PlotOneTrial
    end

    function DrawSpikeWaves(DATA, ispk, nclusters, ctype)
    for j = 1:nclusters+1
        vs{j} = [];
        xs{j} = [];
    end
    splen = size(DATA.AllData.Spikes.values,2);
    adc = DATA.AllData.Spikes.values(ispk,:);
    dvdt = DATA.AllData.Spikes.dVdt(ispk,:);
    for spk = 1:length(ispk);
        j = DATA.AllData.Spikes.codes(ispk(spk), ctype)+1;
        if DATA.plot.dvdt
            vs{j} = [vs{j} dvdt(spk,:) NaN];
            xs{j} = [xs{j} [1:splen-1] NaN];
        else
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
        end
    end
    set(0,'CurrentFigure',DATA.svfig);
    nc = min([nclusters+1 length(DATA.svh)]);
    for j = 1:nc
        if ~isempty(xs{j})
            if ishandle(DATA.svh(j))
                set(DATA.svh(j),'Xdata' , xs{j}, 'Ydata', vs{j});
            else
                DATA.svh(j) = line('Xdata' , xs{j}, 'Ydata', vs{j});
            end
        elseif ishandle(DATA.svh(j)) %no spikes, wipe clean
            set(DATA.svh(j),'Xdata' , 0, 'Ydata', 0);
        end
    end
    drawnow;


function [DATA, ispk] = APlotTrialSpikes(DATA, times, colors, nclusters, classify, varargin)

j = 1;
Trial = [];
while j <= length(varargin)
    if strncmpi(varargin{j},'Trial',4)
        j = j+1;
        Trial = varargin{j};
    end
    j = j+1;
end
if DATA.state.recut
    ctype = 2;
else
    ctype = 1;
end
splen = size(DATA.AllData.Spikes.values,2);
% try calculating energy for all spikes in Expt in one step.
ispk = find(DATA.AllData.Spikes.times > times(1) &...
        DATA.AllData.Spikes.times < times(2));
    if ispk
        [cx, DATA] = GetSpikeVals(DATA,ispk, DATA.plot.clusterX, classify);
        DATA.Spikes.cx(ispk) = cx;
  %      [cy, DATA] = GetSpikeVals(DATA,ispk, SPKVARE, classify);
        [cy, DATA] = GetSpikeVals(DATA,ispk, DATA.plot.clusterY, classify);
        DATA.Spikes.cy(ispk) = cy;
        DATA.currentspike = ispk(1);
            adc = DATA.AllData.Spikes.values(ispk,:);
            dvdt = DATA.AllData.Spikes.dVdt(ispk,:);
% recut == 2 means that clusters are not set here, but clusters have
% been defined (previous list), so use those properites.
         if DATA.state.recut == 2
             DATA = SetSpkCodes(DATA,ispk,0);
         end

         
        for j = 1:nclusters+1
            vs{j} = [];
            xs{j} = [];
        end
        for spk = 1:length(ispk);
            j = DATA.AllData.Spikes.codes(ispk(spk), ctype)+1;
            if DATA.plot.dvdt ==2
                vs{j} = [vs{j} dvdt(spk,:) NaN];
                xs{j} = [xs{j} adc(spk,1:splen-1) NaN];
            elseif DATA.plot.dvdt
                vs{j} = [vs{j} dvdt(spk,:) NaN];
                xs{j} = [xs{j} [1:splen-1] NaN];
            elseif DATA.plot.nodc
            vs{j} = [vs{j} adc(spk,:)-mean(adc(spk,:)) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
            else
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
            end
        end
        set(0,'CurrentFigure',DATA.svfig);
        nc = min([nclusters+1 length(DATA.svh)]);
        for j = 1:nc
            if ~isempty(xs{j})
                if ishandle(DATA.svh(j))
                    set(DATA.svh(j),'Xdata' , xs{j}, 'Ydata', vs{j});
                else
                    DATA.svh(j) = line('Xdata' , xs{j}, 'Ydata', vs{j});
                end
            elseif ishandle(DATA.svh(j)) %no spikes, wipe clean
                    set(DATA.svh(j),'Xdata' , 0, 'Ydata', 0);
            end
        end
        if ~isempty(Trial)
            nspk = sum(DATA.AllData.Spikes.codes(ispk,2) == DATA.currentcluster);
            title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f %d/%d spks',abs(Trial.Trial),...
                Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed,nspk,length(ispk)));
        end
        drawnow;

        set(0,'CurrentFigure',DATA.xyfig);
        for j = 0:nclusters
            sp = find(DATA.AllData.Spikes.codes(ispk, ctype) == j);
            plot(cx(sp),cy(sp),...
                '.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
            hold on; %% need this when called from PlotOneTrial
        end
    end
    DATA.minplottime = 0.00;
    if DATA.minplottime > 0.001
        while toc < DATA.minplottime
        end
    end

function PlotLFP(state, Trial, rate)
    if state.uselfp & isfield(Trial,'LFP')
        set(0,'CurrentFigure',state.lfig);
        times = ([1:length(Trial.LFP)]-Trial.lfpo).*rate;
        plot(times,Trial.LFP);
    end

        
                
function spikelist = WhichClusters(top,varargin)

spikelist = [];

if get(findobj(top,'Tag','UseCluster0'),'value')
    spikelist = [spikelist 0];
end
if get(findobj(top,'Tag','UseCluster1'),'value')
    spikelist = [spikelist 1];
end
if get(findobj(top,'Tag','UseCluster2'),'value')
    spikelist = [spikelist 2];
end
if get(findobj(top,'Tag','UseCluster3'),'value')
    spikelist = [spikelist 3];
end
if get(findobj(top,'Tag','UseCluster4'),'value')
    spikelist = [spikelist 4];
end
if isempty(spikelist)
    spikelist = 0;
end


function DATA = CountSpikes(DATA, expid, varargin)

replot = 0;
j = 1;
while j <= nargin -2
    if strncmp(varargin{j},'replot',3)
        replot = 1;
    end
    j = j+1;
end

nt = 1;
if DATA.state.recut
    ctype = 2;
else
    ctype = 1;
end
if  DATA.state.online == 2
return;
end
if DATA.Expts{expid}.gui.classified == 0
   DATA = SetExptSpikes(DATA, expid, 0);
end
Spks = DATA.AllData.Spikes;


spikelist = WhichClusters(DATA.toplevel);

for trial = [DATA.Expts{expid}.Trials]
    ispk = find(Spks.times > trial.Start(1) & Spks.times < trial.End(end)+500);
    ispks = find(ismember(Spks.codes(ispk,ctype),spikelist));
    DATA.Expts{expid}.Trials(nt).Spikes = round(Spks.times(ispk(ispks)) - trial.Start(1));
    spks = DATA.Expts{expid}.Trials(nt).Spikes;
    DATA.Expts{expid}.Trials(nt).count = sum(spks > 500);
    ispks = find(~ismember(Spks.codes(ispk,ctype),spikelist));
    DATA.Expts{expid}.Trials(nt).OSpikes = round(Spks.times(ispk(ispks)) - trial.Start(1));
    if DATA.state.recut
        DATA.Expts{expid}.Trials(nt).Ocodes = Spks.codes(ispk(ispks),2);
    else
        DATA.Expts{expid}.Trials(nt).Ocodes = Spks.codes(ispk(ispks),1);
    end
    nt = nt+1;
end
if replot
    GetFigure(DATA.tag.dataplot);
    args = PlotArgs(DATA, DATA.Expts{expid});
    PlotExpt(DATA.Expts{expid},args{:});
    t = get(get(gca,'title'),'String');
    title([t 'Cl' sprintf(' %d',spikelist)]);
end
DATA.Expts{expid}.gui.counted = sum(spikelist);
DATA.spikelist = spikelist;

function nc = CountClusters(Cluster)
    nc = 0;
if isempty(Cluster)
    nc = 0;
else
    for j = 1:length(Cluster)
        if ~isempty(Cluster{j});
            nc = nc+1;
        end
    end
end

function outname = CombinedName(DATA,eid,cluster)
suff = '';
it = strmatch(DATA.exptypelist{eid(1)},DATA.expstrs,'exact');
if ~isempty(it)
    expname = DATA.expnames{it};
else
    expname = DATA.exptypelist{eid(1)};
end
stimname = strrep(DATA.explist{eid(1)},['.' DATA.exptypelist{eid(1)}],'');
if strfind(stimname,'RC')
    stimname = strrep(stimname,'RC','');
    suff = 'RC';
end
if isempty(strfind(expname,stimname))
    expname = [stimname '.' expname];
end
outname = [strrep(DATA.datafilename,'.mat','.c') num2str(cluster) '.' expname suff '.mat'];

function SetExptClusters(caller,b)
%DATA = combine('getstate');
DATA = GetDataFromFig(caller);
DATA.Expts{DATA.currentexpt}.Cluster = DATA.cluster;
DATA.Expts{DATA.currentexpt}.gui.clustertype = 1;
DATA.Expts{DATA.currentexpt}.gui.classified = 1;
DATA.Expts{DATA.currentexpt}.gui.ncluster = CountClusters(DATA.cluster);
spkvarnames= DATA.spkvarnames;
for j = 1:length(DATA.Expts)
    if isfield(DATA.Expts{j},'Cluster')
        AllClusters{j}.Cluster = DATA.Expts{j}.Cluster;
        AllClusters{j}.ids = [DATA.Expts{j}.Trials(1).id DATA.Expts{j}.Trials(end).id];
        clustertypes(j) = DATA.Expts{j}.gui.clustertype;
    end
    excludelist{j}  = find([DATA.Expts{j}.Trials.Trial] < 0);
end
clid = DATA.AllData.Spikes.codes(:,2);
save(ClusterFile(DATA),'AllClusters','clustertypes','excludelist','clid','spkvarnames');
% now re-do list of su-expts to reflect cut clusters
eid = get(DATA.clst,'value');
DATA = ListSubExpts(DATA,eid,'relist');
set(DATA.toplevel,'UserData',DATA);
cid = findobj('Tag','ClusterIsSet');
set(cid,'value',1);


function DelClusterButton(caller,b)
it = findobj('Tag','Clusterid');
c = get(it,'value');
DeleteCluster(c, caller);

function ClrSpkWin(caller,b)
%DATA = combine('getstate');
if isstruct(caller)
    DATA = caller;
else
    DATA = GetDataFromFig(caller);
end
GetFigure(DATA.tag.clusterxy);
ym = get(gca,'ylim');
xm = get(gca,'xlim');
hold off;
plot(0,0,'+');
set(gca,'ylim',ym);
set(gca,'xlim',xm);
hold on;
DrawClusters(DATA,0);


function cfile = ClusterFile(DATA,varargin)
getonline = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'getonline',4) %Read Clusters from Online mat file
        getonline = 1;
    end
    j = j+1;
end
if DATA.state.online == 0 %%not online data
    if getonline
        [a,b] = splitpath(DATA.datafilename);
        cfile = [b '/OnlineClusters.mat'];
    else
    cfile = strrep(DATA.datafilename,'.mat','.cl.mat');
    end
else
    cfile = [DATA.datafilename '/Clusters.mat'];
end

function cfile = CombinerLst(DATA)
if isfield(DATA,'datafilename') %%not online data
    cfile = strrep(DATA.datafilename,'.mat','.combine.mat');
else
    cfile = [DATA.datafilename '/combine.mat'];
end

function DATA = CalcClusterVars(DATA, ispk)
SPKENERGY=1;
SPKVARE = 2;
if ispk
    adc = DATA.AllData.Spikes.values(ispk,:);
    if isfield(DATA.AllData.Spikes,'dVdt')
        energy  = sum(DATA.AllData.Spikes.dVdt(ispk,:)'.^2);
    else
    energy  = sum(diff(adc').^2);
    end
    svar = var(adc');
    DATA.Spikes.energy(ispk)= energy;
    if DATA.plot.clusterX == SPKENERGY
        DATA.Spikes.cx(ispk)= energy;
    else
        DATA.Spikes.cx(ispk)= GetSpikeVals(DATA, ispk, DATA.plot.clusterX, 1);
    end
    DATA.Spikes.vw(ispk) = svar./energy;
    if DATA.plot.clusterY == SPKVARE
        DATA.Spikes.cy(ispk)= DATA.Spikes.vw(ispk);
    else
        DATA.Spikes.cy(ispk)= GetSpikeVals(DATA, ispk, DATA.plot.clusterY, 1);
    end
end

function [x,y, DATA] = GetClusterSpace(DATA, Expt)


            x = DATA.plot.clusterX;
        y = DATA.plot.clusterY;
    if isfield(Expt,'Cluster') & ~isempty(Expt.Cluster)
        j = 1;
        while j < length(Expt.Cluster) && isempty(Expt.Cluster{j})
            j = j+1;
        end
        if j <= length(Expt.Cluster) & ~isempty(Expt.Cluster{j})
        x = Expt.Cluster{j}.params(1);
       y = Expt.Cluster{j}.params(2);
       DATA.clusterArange = Expt.Cluster{j}.Arange;
       DATA.clusterBrange = Expt.Cluster{j}.Brange;
       if isfield(Expt.Cluster{j},'Erange')
           DATA.clusterErange = Expt.Cluster{j}.Erange;
       end
        end
    else
        x = DATA.plot.clusterX;
        y = DATA.plot.clusterY;
    end
       

function DATA = SetExptSpikes(DATA, expid, show)

times(1) = DATA.Expts{expid}.Trials(1).Start(1);
times(2) = DATA.Expts{expid}.Trials(end).End(end);
ispk = find(DATA.AllData.Spikes.times > times(1) &...
    DATA.AllData.Spikes.times < times(2));
if ~isfield(DATA,'spklist') | isempty(DATA.spklist)
    DATA.spklist = ispk;
end
DATA.Expts{expid}.gui.spkrange = [min(ispk) max(ispk)];

DATA = CalcClusterVars(DATA, ispk);
% if a cluster is set for this expt, use it.
% otherwise use the current one
if isfield(DATA.Expts{expid},'Cluster') & ~isempty(DATA.Expts{expid}.Cluster)
    DATA.cluster = DATA.Expts{expid}.Cluster;
    DATA.Expts{expid}.gui.classified = 1;
elseif isfield(DATA.cluster{1},'x') %% actually is a cluster defined
    DATA.Expts{expid}.gui.classified = 2;
end
if DATA.state.recut
    DATA = SetSpkCodes(DATA,ispk,show);
end


function DATA = DrawClusters(DATA, setfig)
if setfig
    set(0,'CurrentFigure',DATA.xyfig);
end
for j = 1:length(DATA.cluster)
    C = DATA.cluster{j};
    while ~isempty(C)
    if isfield(C,'params')
        if C.params(1) == DATA.plot.clusterX & ...
                C.params(2) == DATA.plot.clusterY
            h = DrawCluster(C,DATA.spkcolor{j+1});
            if(h)
                DATA.cluster{j}.h = h;
            end
            hold on;
        end
    end
    if isfield(C,'Cluster')
        C = C.Cluster;
    else
        C ={};
    end
    end
end

function h =  DrawCluster(cluster, color)

h =0;
if isempty(cluster) | ~isfield(cluster,'x')
    return;
end
if isfield(cluster,'h') & ishandle(cluster.h)
    delete(cluster.h);
end
tmp.r = [cluster.x(2) cluster.y(2)];
tmp.c = [cluster.x(1) cluster.y(1)];
tmp.xrange = cluster.x(3);
tmp.yrange = cluster.y(3);
tmp.angle = -cluster.angle;
tmp.color = color;
tmp = myellipse(tmp);
h = tmp.lasth;

function DATA = PlaySpikes(DATA, expid)

global mousept;

mode = 2;
cw = DATA.plot.cw;
ch = DATA.plot.ch;
rh = ch+10;
[xyfig, isnew] = GetFigure('ClusterPlot');
DATA.xyfig = xyfig;
if isnew
        if ~isfield(DATA,'figpos') | isempty(DATA.figpos{2})
        bp = get(DATA.toplevel,'Position');
        fp = get(xyfig,'Position');
       DATA.figpos{2} =  [bp(1)+bp(3) bp(2) fp(3) fp(4)];
        set(xyfig,'Position',DATA.figpos{2});
        end
    bp = [5 5 40 20];
    cp = bp;
    uicontrol(xyfig,'style','pop','string','1|2|3|4|5|6|7','Position',bp,'Tag','Clusterid');
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', @DensityPlot, ...
'String', 'Dens','Tag','Density', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', @NextList, ...
'String', 'Next', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'style','pushbutton','string','spool','Position',cp,'Tag','SpoolSpikes',...
        'Callback', @SpoolSpikes);

    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Set','Position',bp,'Callback', @SetExptClusters);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw;
    uicontrol(xyfig,'style','CheckBox','Position',bp,'Tag','ClusterIsSet');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*3;
    uicontrol(xyfig,'style','pushbutton','string','Del','Position',bp,'Callback', @DelClusterButton);
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Clr','Position',bp,'Callback', @ClrSpkWin);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*5;
    uicontrol(xyfig,'style','text','string','Max: X','Position',bp);
    bp(1) = bp(1) + bp(3);
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterXrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterXmax');

    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Y','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterYrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterYmax');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Z','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterZrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterZmax');
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*5;
    uicontrol(xyfig,'style','CheckBox','string','auto','Position',bp,'Tag','AutoScale',...
        'value',DATA.plot.autoscale,'Callback',@Update);
    bp(1) = bp(1) + bp(3) + 10;
    tmpdat.parentfigtag = DATA.tag.top;
    set(xyfig,'UserData',tmpdat);
end
hold off;
%if isempty(findobj('Tag',DATA.tag.spikev
 %   sfig = figure('Renderer','painters','Tag','DATA.tag.spikev');
[sfig, isnew] = GetFigure(DATA.tag.spikev);
if isnew
   if ~isfield(DATA,'figpos') | length(DATA.figpos) < 3 | isempty(DATA.figpos{3})
        bp = get(DATA.toplevel,'Position');
        fp = get(sfig,'Position');
       DATA.figpos{3} =  [bp(1) bp(2)-bp(4) fp(3) fp(4)];
        set(sfig,'Position',DATA.figpos{3});
    end

    x = 10;
    c = 5;
    bp = [10 10 40 20];
    uicontrol(sfig,'style','pushbutton','string','>>','Position',bp,'Tag','NextTrial',...
        'Callback', @PlayNextTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','<<','Position',bp,'Tag','LastTrial',...
        'Callback', @PlayLastTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','X','Position',bp,'Tag','CutTrial',...
        'Callback', @CutTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pop','string','1|2|3','Position',bp,'Tag','ChooseTrial',...
        'Callback', @SelectTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','+','Position',bp,'Tag','NextSpike',...
        'Callback', @PlayNextSpike);
    bp(1) = bp(1)+bp(3);
    bp(3)=c*6;
    uicontrol(sfig,'style','pushbutton','string','spool','Position',bp,'Tag','SpoolSpikes',...
        'Callback', @SpoolSpikes);
    bp(1) = bp(1) + bp(3);
    
    bp(3)=c*9;
    uicontrol(gcf,'Style', 'checkbox',...
    'String', 'dVdt', 'Tag', 'dVdt', 'Position', bp,'value',DATA.plot.dvdt,...
    'Callback',@Update);
    tmpdat.parentfigtag = DATA.tag.top;
    set(sfig,'UserData',tmpdat);

    DATA.hline = 0;    
end
    it = findobj(sfig,'Tag','ChooseTrial');
    set(it,'string',sprintf('%d|',[DATA.Expts{expid}.Trials.Trial]),'value',1);
    DATA.svfig = sfig;
if DATA.spooling ~= 2 %2 == spool from current trial
    DATA.currenttrial = 0;
end
DATA.ISIpair = 0;
DATA.currentexpt = expid;
        ClrSpkWin(DATA);
if DATA.test.fastplot
    set(DATA.svfig,'renderer','painters');
    ax = findobj(DATA.svfig,'type','ax');
    set(ax,'DrawMode','fast','XTickMode','manual','YTickMode','manual');
    set(ax,'Xlim',[1 46],'Ylim',[-5 5]);
end

if DATA.state.uselfp
    DATA.state.lfig = GetFigure('LFP');
    hold off;
    DATA.Expts{expid} = LoadSpike2LFP(DATA.Expts{expid});
end
colors = mycolors;
Spks = DATA.AllData.Spikes;
tic;
nclusters = 1+max(Spks.codes(:,1));
set(0,'CurrentFigure',xyfig);

reclassify = 0;
cid = findobj('Tag','ClusterIsSet');
if DATA.spooling
    DATA.state.recut = 1;
elseif isfield(DATA.Expts{expid},'Cluster') & DATA.state.recut
    DATA.cluster = DATA.Expts{expid}.Cluster;
    DATA.state.recut = 1;
    if DATA.Expts{expid}.gui.clustertype == 2 %% THe online cut
        set(cid,'value',0);
    else
        set(cid,'value',1);
    end
    if ~isfield(DATA.Expts{expid}.gui,'classified') | DATA.Expts{expid}.gui.classified ~= 1
     DATA = SetExptSpikes(DATA,expid,0);
    end
elseif ~isfield(DATA,'cluster')
    DATA.cluster = {};
elseif DATA.state.recut %No cluster yet defined
   DATA.state.recut = 2;
% when inheriting a cluster, don't inherit any spike ranges
  nclusters = 1+max(Spks.codes(:,2));
  for j = 1:nclusters
      if j <= length(DATA.cluster) & ~isempty(DATA.cluster{j})
      if isfield(DATA.Expts{expid}.gui,'spkrange') & ~isempty(DATA.Expts{expid}.gui.spkrange)
      DATA.cluster{j}.firstspk = DATA.Expts{expid}.gui.spkrange(1);
      DATA.cluster{j}.lastspk = DATA.Expts{expid}.gui.spkrange(2);
      else
      DATA.cluster{j}.firstspk = 0;
      DATA.cluster{j}.lastspk = 0;
      end
      end
  end
    set(cid,'value',0);
end

DATA = DrawClusters(DATA,0);

if DATA.plot.autoscale == 0
  set(gca,'Xlim',DATA.plot.clusterXrange,'Ylim',DATA.plot.clusterYrange); 
end 

set(0,'CurrentFigure',sfig);
hold off;

for j = 1:nclusters+1
 DATA.svh(j) = plot([1:46],[1:46] * 0,'color',DATA.spkcolor{j});
 hold on;
end
if ~isfield(DATA.Expts{expid}.gui,'spks') %haven't loaded this expt yet
    times = [DATA.Expts{expid}.Trials(1).Start(1) DATA.Expts{expid}.Trials(end).End(end)];
    nspk = sum(DATA.AllData.Spikes.times > times(1) & DATA.AllData.Spikes.times < times(2)); 
    if nspk < 10000
        DATA.ptsize = 6;
    elseif nspk < 2000
        DATA.ptsize = 10;
    else
        DATA.ptsize = 4;
    end
elseif length((DATA.Expts{expid}.gui.spks)) > 10000
    DATA.ptsize = 4;
elseif length((DATA.Expts{expid}.gui.spks)) > 1000
    DATA.ptsize = 6;
else
    DATA.ptsize = 10;
end
if DATA.plot.setptsize
        DATA.ptsize = DATA.plot.setptsize;
end

if DATA.plot.dvdt == 2
    set(gca,'Xlim',[-5 5],'Ylim',[-5 5]);
else
    set(gca,'Xlim',[1 46],'Ylim',[-5 5]);
end
DATA.nclusters = nclusters;
nt = 1;
allspks = [];
firstspk = 0;
hline = 0;
if ~isfield(DATA,hline)
    DATA.hline = 0;
end

if length(DATA.cluster) < DATA.currentcluster
    DATA.currentcluster = 1;
end
start = max([DATA.currenttrial 1]);
last = length(DATA.Expts{expid}.Trials);
% negative Trial numbers indicate manual exclusion
uset = find([DATA.Expts{expid}.Trials(start:last).Trial] > 0);
for trial = [DATA.Expts{expid}.Trials(uset).Trial]
    if DATA.state.uselfp
    PlotLFP(DATA.state,DATA.Expts{expid}.Trials(uset(nt)),DATA.Expts{expid}.Header.LFPsamplerate);
    end
    set(0,'CurrentFigure',sfig);
    hold off;
    itrial = find(DATA.AllData.Trialids == trial);
    if mode == 1
        DATA = PlotTrialSpikes(DATA,itrial,colors, clusters);
    elseif mode == 2
        times(1) = DATA.Expts{expid}.Trials(uset(nt)).Start(1);
        times(2) = DATA.Expts{expid}.Trials(uset(nt)).End(end);
        Trial = DATA.Expts{expid}.Trials(uset(nt));
        Trial.ed = GetEval(DATA.Expts{expid},'ed',DATA.currenttrial);
        [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial);
        allspks = [allspks spks'];
        nt = nt+1;
    end
    hold on;
end

if DATA.state.uselfp  % look for calibration spikes

 GetFigure('LFP');
 hold off;
 CalcLFPPulse(DATA.Expts{expid},DATA.AllData,'plot');
 GetFigure('SpikeV');
end
if isempty(allspks)
    DATA.spkrange(1) = 1;
    DATA.spkrange(2) = 1;
else
    DATA.spkrange(1) = min(allspks);
    DATA.spkrange(2) = max(allspks);
    DATA.spklist = allspks;
end
DATA.Expts{expid}.gui.spkrange = DATA.spkrange;
DATA.Expts{expid}.gui.spks = allspks;
DATA.Expts{expid}.gui.s2clusters = 1+max(DATA.AllData.Spikes.codes(allspks,1));
DATA.s2clusters = DATA.Expts{expid}.gui.s2clusters;
FinishXYPlot(DATA);
set(xyfig, 'KeyPressFcn',@KeyPressed);
set(xyfig, 'WindowButtonDownFcn',@ButtonPressed);
set(xyfig, 'WindowButtonMotionFcn',@ButtonDragged);
set(xyfig, 'WindowButtonUpFcn',@ButtonReleased);
axis('manual');
ClearMouse;
DATA.densityplot = 0;
xrange = get(gca,'Xlim');
yrange = get(gca,'Ylim');

if DATA.test.fastplot
    set(DATA.xyfig,'renderer','painters');
    ax = findobj(DATA.xyfig,'type','ax');
    set(ax,'DrawMode','fast','XTickMode','manual','YTickMode','manual');
    %set(,'DrawMode','fast');
end
if DATA.state.fixrange
    if length(allspks) > 1000
        prc = 99;
    else
        prc = 95;
    end
    xm = prctile(DATA.Spikes.energy(allspks),prc);
    ym = prctile(DATA.Spikes.vw(allspks),prc);
    if yrange(2) > 2 * ym
        set(gca,'Ylim',[0 ym * 1.5]);
    end
    if xrange(2) > 2 * xm
        set(gca,'Xlim',[0 xm * 1.5]);
    end
end
mousept.xrange = diff(xrange);
mousept.yrange = diff(yrange);
toc;
SetGui(DATA);

function DATA = GetDataFromFig(a)

    DATA = get(a,'UserData');

    if isempty(DATA)
        b = get(a,'parent');
        DATA = get(b,'UserData');
    end
    if isfield(DATA,'parentfigtag')
        DATA = get(findobj('Tag',DATA.parentfigtag),'UserData');
    end
    
function RescaleClusterPlot(a,b)
    DATA = GetDataFromFig(a);
it = findobj(DATA.xyfig,'Tag','ClusterXmax');
if it
    DATA.plot.clusterXrange(2) = str2num(get(it,'string'));
end

it = findobj(DATA.xyfig,'Tag','ClusterYmax');
if it
    DATA.plot.clusterYrange(2) = str2num(get(it,'string'));
end
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
if it
    DATA.plot.clusterZrange(2) = str2num(get(it,'string'));
end
ax = findobj(DATA.xyfig,'type','ax');

%
%if the user has typed in a value, surely they want to manually scales
DATA.plot.autoscale = 0;
set(findobj(DATA.xyfig,'Tag','AutoScale'),'value',DATA.plot.autoscale);
if DATA.plot.autoscale == 0
    set(ax,'Xlim',DATA.plot.clusterXrange);
    set(ax,'Ylim',DATA.plot.clusterYrange);
end
if DATA.densityplot
    caxis([0 DATA.plot.clusterZrange(2)]);
end
set(DATA.toplevel,'UserData',DATA);

function SpoolSpikes(a,b)
DATA = GetDataFromFig(a);
DATA.spooling = 1;
if strmatch(get(a,'Tag'),'SpoolSpikes')
    DATA.spooling = 2; %spool from current spike to end
end
DATA = PlaySpikes(DATA,DATA.currentexpt);
set(DATA.toplevel,'UserData',DATA);

function NextList(a,b)
%DATA = combine('getstate');
DATA = GetDataFromFig(a);

strs = get(DATA.elst,'string');
val = get(DATA.elst,'value');
if val < length(strs)
    set(DATA.elst,'value',val+1)
    combine('setexpt',DATA);
end


function ClearMouse()
global mousept;

mousept.start = [];
mousept.lasth = [];
mousept.angle = 0;
mousept.r = [1 1];
mousept.down = 0;
mousept.mode = 0;


function DensityPlot(a,b)
global mousept;

%DATA = combine('getstate');
DATA = GetDataFromFig(a);

if isfield(DATA,'spklist') & length(DATA.spklist) > 10
    expspks = DATA.spklist;
else
    expspks = DATA.spkrange(1):DATA.spkrange(2);
end

    
    figure(DATA.xyfig);
if DATA.densityplot
    DATA.densityplot = 0;
    hold off;
    DATA = DrawXYPlot(DATA,expspks);
%    ClearMouse;
    set(DATA.toplevel,'UserData',DATA);
    SetGui(DATA);
    return;
end

energy = DATA.Spikes.cx(expspks);
vw = DATA.Spikes.cy(expspks);

if length(vw) > 10000
    lprc = 0.01;
    hprc = 99.99;
elseif length(vw) > 10000
    lprc = 0.1;
    hprc = 99.9;
    
elseif length(vw) > 1000
    lprc = 1;
    hprc = 99;
else
    lprc = 5;
    hprc = 95;
end    

erange = [prctile(energy,lprc) prctile(energy,hprc)];
vrange = [prctile(vw,lprc) prctile(vw,hprc)];
erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
%GetFigure('DensityPlot');
hold off;
nbins = 200;
[x,y] = meshgrid(linspace(erange(1),erange(2),nbins),linspace(vrange(1),vrange(2),nbins));
mode = 2;
tic;
if mode ==1 % add real gaussian to grid for each
    z = zeros(size(x));
    sx = (diff(erange)/100)^2;
    sy = (diff(vrange)/100)^2;
    for j=1:length(energy)
        z = z + exp(-(x-energy(j)).^2/sx - (y-vw(j)).^2/sy);
    end
elseif mode ==2 %build fine 2-D histogram, then smooth
    [gx,gy] = meshgrid(-10:10,-10:10);
    sx=3;
    sy=3;
    G = exp(-(gx).^2/sx - (gy).^2/sy);
    G = G./sum(G(:));
    z = zeros(size(x));
    vi = 1+floor(nbins * (vw-vrange(1))/diff(vrange));
    ei = 1+floor(nbins * (energy-erange(1))/diff(erange));
    idx = find(ei > 0 & ei <= nbins & vi > 0 & vi <= nbins);
    for j =idx
        z(vi(j),ei(j)) = z(vi(j),ei(j))+1;
    end
    z = conv2(z,G,'same');
end
toc
pcolor(x,y,z);
shading('interp')
hold on;
DATA = DrawClusters(DATA, 0);
DATA.densityplot = 1;
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
x = caxis;
set(it,'string',sprintf('%.2f',x(2)));
SetGui(DATA);
set(DATA.toplevel,'UserData',DATA);


function DATA = ListExpts(DATA, Expts);

SpkDefs;
explist = {};
na = 1;
nb = 1;

explist{1} = 'All';
explabel{1} = 'All';

exptypelist = [];
na= 2;
for j = 1:length(Expts)
    stimname = stimnames{Expts{j}.Stimvals.st+1};
    if strmatch(Expts{j}.Stimvals.e2, 'e0')
        if strmatch(Expts{j}.Stimvals.e3, 'e0')
        exptypename = Expts{j}.Stimvals.et;
        else
            exptypename = [Expts{j}.Stimvals.et 'X' Expts{j}.Stimvals.e3];
        end
            
%        expname = [Expts{j}.Stimvals.et];
    else
        exptypename = [Expts{j}.Stimvals.et 'X' Expts{j}.Stimvals.e2];
    end
    if isfield(Expts{j}.Header,'Options') & strfind(Expts{j}.Header.Options,'+cr')
        exptypename = ['C' exptypename];
    end
    if Expts{j}.Header.rc
 %       exptypename = [exptypename 'RC'];
    end
    expname = [stimname '.' exptypename];
% if this expt is not in the list, add it
    if isempty(strmatch(Expts{j}.Header.expname,explist,'exact'))
        explist{na} = expname;
        if DATA.state.online == 2
            explist{na} = Expts{j}.Header.expname;
%            Expts{j}.Header.expname = expname;
        else
            explist{na} = Expts{j}.Header.expname;
        end
        exptypelist{na} = exptypename;
        DATA.explist = explist;
        DATA.exptypelist = exptypelist;
        cls = '';
        for k = 1:3
            outname = CombinedName(DATA, na,k);
            if exist(outname,'file')
                cls = [cls ' c' num2str(k)];
            end
        end
        explabel{na} = [explist{na} cls];
        DATA.explist = explist;
        na = na+1;
    else
%        explist{na} = 'unknown';
    end    
end
p  = get(DATA.clst,'Listboxtop');
if p > length(explabel)
    set(DATA.clst,'ListboxTop',1);
end
p  = get(DATA.clst,'value');
if p > length(explabel)
    set(DATA.clst,'value',1);
end
set(DATA.clst,'string',explabel);
DATA.explist = explist;
DATA.exptypelist = exptypelist;

function [id, ok] = CheckListForName(list,name)

fid = fopen(list,'r');
if fid >0 
    ok = 1;
    names = textscan(fid,'%s');
    files = splitpath(names{1});
    id = strmatch(splitpath(name),files);
    fclose(fid);
else
    ok = 0;
    id = [];
end


function CheckLists(DATA)
       
 pairs = {'ORBW' '/bgc/bgc/anal/orbw/lemORBW.lst' ; ...
     'OTRC' '/bgc/bgc/anal/orbw/otrc.lst'};
        
        [a,dp] = splitpath(DATA.datafilename);
        d = dir(dp);
        nf = 0;
        for j = 1:length(d)
            if strfind(d(j).name,'.mat') 
              for k = 1:size(pairs,1)
                if strfind(d(j).name,pairs{k,1})
                    if isempty(CheckListForName(pairs{k,2},d(j).name))
                    nf = nf+1;
                    dat.files{nf} = [dp '/' d(j).name];
                    dat.lists{nf} = pairs{k,2};
                    else
                        fprintf('%s Already in %s\n',d(j).name,pairs{k,2});
                    end
                    dat.nf = nf;
                end
              end
            end
        end
if nf
    SPACE = 10;
    cw=8;
    ch=10;
    figure('Position',[10 10 cw*60, (ch+SPACE) * (nf+2)]);
    fp = get(gcf,'Position');
    bp = [10 10 fp(4) 20];
    for j = 1:nf
        bp = [10 bp(2)+bp(4)+SPACE fp(3) ch * 2];
        uicontrol(gcf,'Style', 'checkbox',...
            'String', [dat.files{j} dat.lists{j}], 'Tag', 'SuffList', 'Position', bp,'UserData',j);
    end
    bp = [10 10 cw * 5 20];
        uicontrol(gcf,'Style', 'pushbutton',...
            'String','go', 'callback',@cplists);
     set(gcf,'UserData',dat);
    
end


function cplists(caller,b)

    lfig = get(caller,'parent');
    dat = get(lfig,'UserData');
    it = findobj(lfig,'Tag','SuffList');
    for j = 1:length(it)
        go = get(it(j),'value');
        if go
            fprintf('Adding %s to %s\n',dat.files{j},dat.lists{j});
            AddToList(dat.lists{j},dat.files{j});
        end
    end
    close(lfig);
 
function DATA = ListSubExpts(DATA, id, varargin)

setv = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'relist',4)
        setv = get(DATA.elst,'value');
    end
    j = j+1;
end
suffs = DATA.suffs;
na = 1;
subexplist = {};
expid = [];
nrp = zeros(size(DATA.explist));
for j = 1:length(DATA.Expts);
    eid = strmatch(DATA.Expts{j}.Header.expname,{DATA.explist{id}},'exact');
    if( id == 1 | ~isempty(eid)) & isfield(DATA.Expts{j},'Trials')
        tid = regexp(DATA.Expts{j}.Header.Name,'Expt[0-9]*.mat');
        if ~isempty(tid) %online file
            label = sprintf(' (%s:%d-%d)',DATA.Expts{j}.Header.Name(tid:end-4),...
                DATA.Expts{j}.Trials(1).Trial,...
                DATA.Expts{j}.Trials(end).Trial');
        else
            label = sprintf(' %d-%d (id %d-%d)',DATA.Expts{j}.Trials(1).Trial,...
                DATA.Expts{j}.Trials(end).Trial,DATA.Expts{j}.Trials(1).id,DATA.Expts{j}.Trials(end).id);
        end
        if DATA.Expts{j}.Header.psych
            label = [label ' P'];
        end
        if DATA.Expts{j}.Header.rc
            label = [label ' RC'];
        end
        if id == 1
            expi = strmatch(DATA.Expts{j}.Header.expname, DATA.explist,'exact');
            nrp(expi) = nrp(expi)+1;
            subexplist{na} = [DATA.Expts{j}.Header.expname num2str(na) label];
        else
            subexplist{na} = [DATA.explist{id(eid)} suffs(mod(na-1,length(suffs))+1) label];
            subexplist{na} = [DATA.explist{id(eid)} num2str(na) label];
        end
        DATA.explabels{j} = subexplist{na};
        if DATA.Expts{j}.gui.clustertype == 0
            subexplist{na} = [subexplist{na} '*'];
        elseif DATA.Expts{j}.gui.clustertype == 2
            subexplist{na} = [subexplist{na} '(O)'];
        elseif DATA.Expts{j}.gui.ncluster == 0
            subexplist{na} = [subexplist{na} '(Z)'];            
        end
        expid(na) = j;
        na = na+1;
    elseif strmatch('unknown',{DATA.explist{id}})
        subexplist{na} = [DATA.explist{id} suffs(na)];
        expid(na) = j;
        na = na+1;        
    end
end
p  = get(DATA.elst,'Listboxtop');
if p > length(subexplist)
    set(DATA.elst,'ListboxTop',1);
end
    set(DATA.elst,'string',subexplist,'value',setv);
    DATA.expid = expid;
    DATA.subexplist = subexplist;


    
function timerfna(varargin)
DATA = get(findobj('Tag','Combiner'),'UserData');
if DATA.state.autolist
tic;
    combine('relist');
%toc
end

function timerfn(varargin)
    DATA = get(findobj('Tag','Combiner'),'UserData');
    d = dir(DATA.datafilename);
    if d.datenum > DATA.lastread || d.bytes > DATA.lastsize;
 %       fprintf('Reading from %d\n',DATA.linesread);
        [Expts, DATA.linesread] = ReadOnlineTxt(DATA.datafilename, DATA);
%        fprintf('Read %d expts\n',length(Expts));
         Expts = CountTxtSpikes(Expts,DATA.probe,DATA.spikelist);
%         fprintf('List\n');
         if length(Expts) > length(DATA.Expts)
             DATA.Expts = Expts;
             DATA = ListExpts(DATA,Expts);
 %            fprintf('listexp\n');
             DATA = combine('listexps',DATA,'Tag',DATA.tag.top);
%             fprintf('Set\n');
             set(DATA.elst,'value',length(get(DATA.elst,'string')))
             fprintf('Replot (%d,%d) at %s\n',length(Expts),DATA.linesread,datestr(now));
         else
             DATA.Expts = Expts;
             fprintf('Replot (%d,%d,%d) at %s\n',length(Expts),length(Expts{end}.Trials),DATA.linesread,datestr(now));
         end
        DATA.lastread = d.datenum;
        DATA.lastsize = d.bytes;
        combine('setexp',DATA,'Tag',DATA.tag.top);
%        fprintf('Expt Set\n');
    end
    
    
function DATA = BuildGUI(DATA)

scrsz = get(0,'Screensize');
cw= scrsz(3)/140;
ch= scrsz(4)/60;
if scrsz(3) > 2000
    cw = 10;
end
if scrsz(4) > 1200
    ch = 10;
end
DATA.plot.cw = cw;
DATA.plot.ch = ch;
wsiz = [cw*40,ch*40];
SPACE = ch;

cntrl_box = figure('Position', [100 scrsz(4)-(wsiz(2)+ch*4) wsiz(1) wsiz(2)],...
    'NumberTitle', 'off', 'Tag',DATA.tag.top,'Name',DATA.tag.top);
bp = [10 10 cw*38 ch * 10];
DATA.elst = uicontrol(gcf, 'Style','listbox',...
    'Callback', ['combine(''setexpt'',''Tag'',''' DATA.tag.top ''');'],'Tag','subexptlist',...
    'Position',bp,'Max',3,'Min',1);
bp = [10 bp(2)+bp(4)+SPACE cw*7 ch * 2];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['combine(''combine'',''Tag'',''' DATA.tag.top ''');'],...
'String', 'Combine', 'Position', bp);
bp = [bp(1)+bp(3)+SPACE/2 bp(2) cw*5 ch * 2];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['combine(''Save'',''Tag'',''' DATA.tag.top ''');'],...
'String', 'Save', 'Position', bp);
bp = [bp(1)+bp(3)+SPACE/2 bp(2) cw*30 ch * 2];
DATA.saveitem = uicontrol(gcf,'Style','Edit','String','Save As','Position',bp);

bp = [10 bp(2)+bp(4)+SPACE cw*38 ch * 10];
DATA.clst = uicontrol(gcf, 'Style','listbox',...
    'Callback', ['combine(''listexps'',''Tag'',''' DATA.tag.top ''');'],'Tag','explist',...
    'Position',bp,'Max',3,'Min',1);

bp = [10 bp(2)+bp(4)+SPACE cw*9 ch * 2];
uicontrol(gcf,'Style','Text','String','Data File','Position',bp);
bp = [bp(1)+bp(3)+SPACE bp(2) cw*30 ch * 2];
uicontrol(gcf,'Style','Edit','String',DATA.datafilename,'Position',bp,'Callback',['combine(''newfile''''Tag'',''' DATA.tag.top ''');'],...
    'Tag','FileName');

bp = [10 bp(2)+bp(4)+SPACE cw*9 ch * 2];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'ShowN', 'Tag', 'ShowN', 'Position', bp);

bp(1) = bp(1) + bp(3) + SPACE;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Spikes', 'Tag', 'ShowSpikes', 'Position', bp);
bp(1) = bp(1) + bp(3) + SPACE;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Recut', 'Tag', 'Recut', 'Position', bp,'value',DATA.state.recut,...
'Callback',@Update);

bp(1) = bp(1) + bp(3) + SPACE;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Auto', 'Tag', 'AutoPlot', 'Position', bp,'value',DATA.state.autoplot,...
'Callback',@Update);

bp = [10 bp(2)+bp(4)+SPACE cw*6 ch * 2];
uicontrol(gcf,'Style', 'pop',...
'String', 'Mean|Seq:time|Seq:Trials|Seq:Id|Seq:Edit', 'Tag', 'PlotSeq', 'Position', bp,'value',DATA.state.plotseq+1,...
'Callback',@Update);

bp = [bp(1)+bp(3)+SPACE bp(2) cw*6 ch * 2];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Recount', 'Tag', 'Recount', 'Position', bp,'value',DATA.state.recount,...
'Callback',@Update);
bp = [bp(1)+bp(3)+SPACE bp(2) cw*6 ch * 2];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Psych', 'Tag', 'PlotPsych', 'Position', bp,'value',DATA.state.plotpsych,...
'Callback',@Update);
bp = [bp(1)+bp(3)+SPACE bp(2) cw*6 ch * 2];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Combined', 'Tag', 'PlotCombined', 'Position', bp,'value',DATA.state.plotcombined,...
'Callback',@Update);



bp = [10 bp(2)+bp(4)+SPACE cw*3 ch * 2];
uicontrol(gcf,'Style', 'checkbox',...
'String', '0', 'Tag', 'UseCluster0', 'Position', bp,'Callback',{@SetClusters, DATA.tag.top});

bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '1', 'Tag', 'UseCluster1', 'Position', bp,'Callback',{@SetClusters, DATA.tag.top});

bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '2', 'Tag', 'UseCluster2', 'Position', bp,'Callback',{@SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '3', 'Tag', 'UseCluster3', 'Position', bp,'Callback',{@SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '4', 'Tag', 'UseCluster4', 'Position', bp,'Callback',{@SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @FitButton, ...
'String', 'Fit', 'Tag','FitButton','Position', bp);
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'combine(''setexpplot'')', ...
'String', 'Plot', 'Tag','PlotButton','Position', bp);
bp(1) = bp(1)+bp(3);
bp(3) = cw*5;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'LFP', 'Tag', 'UseLFP', 'Position', bp,'Callback',@Update);

bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'pop','String',DATA.probenames,'Position', bp,...
      'Tag','ProbeId','Callback',@Update,'value',1);

  hm = uimenu(gcf,'Label','Mark');
  uimenu(hm,'Label','&Mark','Callback','combine(''Mark'')');
  uimenu(hm,'Label','&relist','Callback','combine(''relist'')');
  uimenu(hm,'Label','&Options','Callback','combine(''options'')');
  uimenu(hm,'Label','&ShowVals','Callback','combine(''showvals'')');
  uimenu(hm,'Label','&To Front','Callback','combine(''winfront'')');
  uimenu(hm,'Label','&Update RF','Callback','combine(''rfupdate'')');
  uimenu(hm,'Label','&Next CEell','Callback','combine(''nextcell'')');
  uimenu(hm,'Label','&Update Lists','Callback','combine(''checklists'')');
  uimenu(hm,'Label','&Close','Callback',['combine(''Close'',''Tag'',' '' DATA.tag.top '' ')']);

  set(gcf,'Menubar','none');
DATA.toplevel = cntrl_box;

function FitButton(a,b)

DATA = GetDataFromFig(a);
[DATA.Expt, plotres] = PlotCombined(DATA, DATA.Expt);
fit = FitExpt(plotres,'plot');
type = strmatch(plotres.type{1},{'Op' 'Pp'});
if ~isempty(type)
    if type(1) == 1
        DATA.fitvals.Op = fit.mean;
        DATA.fitvals.Opw = fit.sd;
        DATA.fitvals.OpRo = GetEval(DATA.Expt,'Ro');
        DATA.Expt.fits.Op = fit;
    elseif type(1) == 2
        DATA.fitvals.Pp = fit.mean;
        DATA.fitvals.PpRo = GetEval(DATA.Expt,'Ro');
        DATA.fitvals.Ppw = fit.sd;
        DATA.Expt.fits.Pp = fit;
    end
end
set(DATA.toplevel,'UserData',DATA);

function cntrl_box = setshow(DATA, tag)
wsc = DATA.wsc;
SPACE = 3 * wsc;
VSPACE = 5 * wsc;
h = 220 * wsc;
w = 350 * wsc;
scrsz = get(0,'Screensize');
cw = DATA.plot.cw;
ch = DATA.plot.ch;
SpkDefs;
   bp(1) = SPACE;
   bp(2) = ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
cntrl_box = figure('Position', [200 scrsz(4)-(h+30)*wsc w*wsc h*wsc], 'Menubar', 'none',...
    'NumberTitle', 'off', 'Tag',tag,'Name','Showvals');
fn = fields(DATA.show);
for j = 1:length(fn)
    id = strmatch(fn{j},CodeNames.Codes);
    uicontrol(gcf,'Style', 'CheckBox','String',CodeNames.Label{id},'Position', bp,...
   'Tag',fn{j},'Callback',@ShowUpdate,'value',DATA.show.(fn{j}));
bp(2) = bp(2)+ch;
end

function str = vec2str(x)
          str = [int2str(x(1)) ':' int2str(x(end))];

function cntrl_box = setoptions(DATA, tag)

wsc = DATA.wsc;
SPACE = 3 * wsc;
VSPACE = 5 * wsc;
h = 220 * wsc;
w = 350 * wsc;
ch = DATA.plot.ch;
cw = DATA.plot.cw;
scrsz = get(0,'Screensize');

cntrl_box = findobj('Tag',tag,'Name','Options');
if ~isempty(cntrl_box)
    figure(cntrl_box);
    return;
end
if ~isfield(DATA,'figpos') | isempty(DATA.figpos{1})
   bp = get(DATA.toplevel,'Position');    
   DATA.figpos{1} =  [bp(1)+bp(3) bp(2) w*wsc h*wsc]
end
dat.parentfigtag = DATA.tag.top;
cntrl_box = figure('Position', DATA.figpos{1} , 'Menubar', 'none',...
    'NumberTitle', 'off', 'Tag',tag,'Name','Options','UserData',dat);


top = num2str(DATA.toplevel); 
bp(1) = SPACE;
bp(2) = ch+VSPACE;
bp(3) = cw*9;
bp(4) = ch+VSPACE;
   cw = 10 * wsc;
   ch = 11 * wsc;
   bh = 18*wsc;
   bp(1) = SPACE;
   bp(2) = ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
   if ~isfield(DATA.plot,'acov')
       DATA.plot.acov = 0;
       DATA.plot.collapse = 0;
       DATA.plot.flip = 0;
   end
   uicontrol(gcf,'Style', 'CheckBox','String','Limit Range','Position', bp,...
      'Tag','FixRange','Callback',@Update,'value',DATA.state.fixrange);
   bp(2) = bp(2)+ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','ShowEM','Position', bp,...
      'Tag','ShowEM','Callback',@Update,'value',DATA.plot.showem);


   bp(2) = bp(2)+ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','Flip','Position', bp,...
      'Tag','Flip','Callback',@Update,'value',DATA.plot.flip);
   bp(2) = bp(2)+ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','Collapse 1','Position', bp,...
      'Tag','Collapse1','Callback',@Update,'value',DATA.plot.collapse);
   bp(2) = bp(2)+ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','Acov','Position', bp,...
      'Tag','Acov','Callback',@Update,'value',DATA.plot.acov);
   bp(2) = bp(2)+ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','ISIH','Position', bp,...
      'Tag','ISIH','Callback',@Update,'value',DATA.plot.showISI);
   bp(1) = bp(1)+bp(3)+SPACE;
   uicontrol(gcf,'Style', 'pushbutton','String','Plot ISIH','Position', bp,...
      'Callback','combine(''PlotISI'')');

  bpa = bp;
  % second column of checkboxes 
  bh = 18*wsc;
   bp(1) = SPACE+bp(3);
   bp(2) = ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','dV vs V','Position', bp,...
      'Tag','PhasePlot','Callback',@Update,'value',(DATA.plot.dvdt == 2));
   bp(2) = bp(2)+ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','Remove DC','Position', bp,...
      'Tag','RemoveDC','Callback',@Update,'value',DATA.plot.nodc);

   bp(2) = bp(2)+ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','Auto relist','Position', bp,...
      'Tag','AutoList','Callback',@Update,'value',DATA.state.autolist);

  bp(2) = bp(2)+ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
   uicontrol(gcf,'Style', 'CheckBox','String','AutoFit','Position', bp,...
      'Tag','AutoFit','Callback',@Update,'value',DATA.state.autofit);

  bp = bpa;
   
   bp(1) = SPACE;
  bp(2) = bp(2)+ch+VSPACE;
   uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'Position', bp,...
      'Tag','ClusterX','Callback',@Update,'value',DATA.plot.clusterX);
   
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.clusterXrange(1)),'Position', bp,...
      'Tag','ClusterXmin','Callback',@Update);


  bp(2) = bp(2)+ch+VSPACE;
   bp(1) = SPACE;
   bp(3) = cw*9;
   uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'Position', bp,...
      'Tag','ClusterY','Callback',@Update,'value',DATA.plot.clusterY);

  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.clusterYrange(1)),'Position', bp,...
      'Tag','ClusterYmin','Callback',@Update);

  bp(2) = bp(2)+ch+VSPACE;
   bp(1) = SPACE;
  uicontrol(gcf,'Style', 'text','string','Pt Range A','Position',bp);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'edit','string',vec2str(DATA.clusterArange),'Position', bp,...
      'Tag','ClusterArange','Callback',@Update);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw;
  uicontrol(gcf,'Style', 'text','string','B','Position',bp);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'edit','string',vec2str(DATA.clusterBrange),'Position', bp,...
      'Tag','ClusterBrange','Callback',@Update);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw;
  uicontrol(gcf,'Style', 'text','string','E','Position',bp);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'edit','string',vec2str(DATA.clusterErange),'Position', bp,...
      'Tag','ClusterErange','Callback',@Update);

  bp(2) = bp(2)+ch+VSPACE;
   bp(1) = SPACE;
   uicontrol(gcf,'Style', 'text','string','N Min','Position',bp);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.nmin),'Position', bp,...
      'Tag','Nmin','Callback',@Update);

    bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'text','string','CP','Position',bp);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'pop','String','None|CP-time|CP trials|CP-hist','Position', bp,...
      'Tag','ShowCP','Callback',@Update,'value',DATA.plot.showcp+1);

    bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'text','string','PtSz','Position',bp);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5; 
  uicontrol(gcf,'Style', 'pop','String','Auto|1|2|3|4|5|6|7|8', 'Position',bp,...
      'Tag','SetPtSize','Callback',@Update,'value',DATA.plot.setptsize+1);

  
  bp(2) = bp(2)+ch+VSPACE;
   bp(1) = SPACE;
   uicontrol(gcf,'Style', 'text','string','sdfw','Position',bp);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 5;
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.sdfw),'Position', bp,...
      'Tag','Sdfw','Callback',@Update);

      
  function ShowUpdate(a,b)

DATA = GetDataFromFig(a);
fn = fields(DATA.show);
for j = 1:length(fn)
    DATA.show.(fn{j}) = GetShow(DATA,fn{j});
end
set(DATA.toplevel,'UserData',DATA);

function [value, it] = GetShow(DATA,tag)
it = findobj(DATA.showid, 'Tag',tag);
if ~isempty(it) 
    value = get(it(1),'value');
else
    value = 0;
end


  

function SetGui(DATA);
    SetCheck('Recut',DATA.state.recut > 0); %% in case it is 2
    SetCheck('Recount',DATA.state.recount);
    SetCheck('UseCluster1',ismember(1,DATA.spikelist));
    if isfigure(DATA.xyfig)
    it = findobj(DATA.xyfig,'Style', 'pushbutton', 'Tag','Density');
    if it
        if DATA.densityplot
            set(it,'String','Pts');
        else
            set(it,'String','Dens');
        end
    
    ax = findobj(DATA.xyfig,'Type','axes');
    if DATA.plot.autoscale
        set(ax,'Ylimmode','auto','Xlimmode','auto');
        DATA.plot.clusterXrange  = get(ax,'Xlim');
        DATA.plot.clusterYrange  = get(ax,'Ylim');
        SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
        SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
    else
        set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
    end
    end
    end
    if isfigure(DATA.optionfig)
    SetField(DATA.optionfig,'ClusterXmin',DATA.plot.clusterXrange(1));
    SetField(DATA.optionfig,'ClusterYmin',DATA.plot.clusterYrange(1));
    end


function SetClusters(a,b,tag)
DATA = get(findobj('Tag',tag),'UserData');
DATA.spikelist = WhichClusters(DATA.toplevel);
if DATA.state.online == 2
    DATA.Expts = CountTxtSpikes(DATA.Expts,DATA.probe,DATA.spikelist);
end
if DATA.state.autoplot
    combine('setexpt','Tag', tag);
end
set(DATA.toplevel,'UserData',DATA);

function Update(a,b)

DATA = GetDataFromFig(a);
DATA = combine('getstate');

if DATA.xyfig & get(a, 'Parent') == DATA.xyfig
    DATA.plot.autoscale = GetCheck('AutoScale',DATA.xyfig);
    it = findobj(DATA.xyfig,'Tag','ClusterId');
    DATA.currentcluster = get(it,'value');
        ax = findobj(DATA.xyfig,'Type','axes');
    if DATA.plot.autoscale
        set(ax,'Ylimmode','auto','Xlimmode','auto');
        DATA.plot.clusterXrange  = get(ax,'Xlim');
        DATA.plot.clusterYrange  = get(ax,'Ylim');
        SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
        SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
    else
        set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
    end
else
DATA.plot.dvdt = GetCheck('dVdt');
DATA.plot.nodc = GetCheck('RemoveDC');
dvv = GetCheck('PhasePlot');
if DATA.plot.dvdt && dvv
    DATA.plot.dvdt = 2;
end

DATA.state.recut = GetCheck('Recut');
DATA.state.recount = GetCheck('Recount');
DATA.state.plotpsych = GetCheck('PlotPsych');
DATA.state.plotcombined = GetCheck('PlotCombined');
DATA.state.showspikes = GetCheck('ShowSpikes');
[DATA.state.autoplot, h] = GetCheck('AutoPlot');
DATA.spikelist = WhichClusters(DATA.toplevel);
DATA.state.uselfp = get(findobj(DATA.toplevel,'Tag','UseLFP'),'value');
id = get(findobj(DATA.toplevel,'Tag','ProbeId'),'value');
DATA.probe= DATA.probelist(id);
DATA.state.plotseq = get(findobj(DATA.toplevel,'Tag','PlotSeq'),'value')-1;


if DATA.state.online %%no LFP available
    DATA.state.uselfp = 0;
end
id = regexp(DATA.outname,'.c[0-9].');
if id
    DATA.outname(id+2) = num2str(DATA.spikelist(1));
    set(DATA.saveitem,'string',DATA.outname);
end
it = findobj('Tag',DATA.tag.options);
if ~isempty(it)
    DATA.state.fixrange = GetCheck('FixRange');
    DATA.plot.showem = GetCheck('ShowEM',it);
    DATA.plot.showcp = get(findobj(it,'Tag','ShowCP'),'value')-1;
    DATA.plot.clusterX = get(findobj(it,'Tag','ClusterX'),'value');
    DATA.plot.clusterY = get(findobj(it,'Tag','ClusterY'),'value');
    DATA.plot.clusterXrange(1) = GetField('ClusterXmin',it);
    DATA.plot.clusterYrange(1) = GetField('ClusterYmin',it);
    DATA.clusterArange = GetField('ClusterArange',it);
    DATA.clusterBrange = GetField('ClusterBrange',it);
    DATA.clusterErange = GetField('ClusterErange',it);
    DATA.plot.nmin = GetField('Nmin',it);
    DATA.plot.sdfw = GetField('Sdfw',it);
    if length(DATA.plot.clusterX) > 1 || length(DATA.plot.clusterY) > 1
        fprintf('ClusterX/Y is too big');
    end
    DATA.plot.acov = GetCheck('Acov',it);
    DATA.plot.flip = GetCheck('Flip',it);
    DATA.plot.collapse = GetCheck('Collapse1',it);
    DATA.plot.showISI = GetCheck('ISIH',it);
    [DATA.state.autolist, h] = GetCheck('AutoList');
    DATA.plot.setptsize = get(findobj(it,'Tag','SetPtSize'),'value')-1;
end
end
set(DATA.toplevel,'UserData',DATA);

if strmatch(get(a,'Tag'),{'PlotSeq' 'Psych'})
    PlotCombined(DATA,DATA.Expt);
end
    if strmatch(get(a,'Tag'),{'ClusterX' 'ClusterY'})
        if isfield(DATA,'spklist') && ~isempty(DATA.spklist)
        expspks = DATA.spklist;
    else
        expspks = DATA.spkrange(1):DATA.spkrange(2);
        end
    GetFigure(DATA.xyfig);
    hold off;
    DATA = CalcClusterVars(DATA, expspks);
    if DATA.plot.setptsize
        DATA.ptsize = DATA.plot.setptsize;
    end
    DrawXYPlot(DATA,expspks);
end

if DATA.state.autoplot & a ~= h & ~ DATA.state.showspikes
    combine('setexp');
end
    
function SetField(parent, tag, value, varargin)
if isfigure(parent)
    it = findobj(parent,'Tag', tag);
else
    it = findobj('Tag', tag);
end
if it
    set(it,'string',num2str(value));
end

function value = GetField(tag, varargin)

if nargin == 2 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    value = str2num(get(it(1),'string'));
else
    value = NaN;
end

function [value, it] = GetCheck(tag, varargin)

if nargin == 2 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    value = get(it(1),'value');
else
    value = 0;
end

function [success] = SetCheck(tag, value,  varargin)

if nargin == 3 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    set(it(1),'value',value);
    success = 1;
else
    success = 0;
end
    
function args = PlotArgs(DATA, Expt)

args = {};
psych = get(findobj('Tag','PlotPsych','Parent',DATA.toplevel),'value');
seq = get(findobj('Tag','PlotSeq','Parent',DATA.toplevel),'value')-1;
on = get(findobj('Tag','ShowN','Parent',DATA.toplevel),'value');
if on
    args = {args{:} 'shown'};
end
on = get(findobj('Tag','PlotSeq','Parent',DATA.toplevel),'value')-1;
if psych & seq
        args = {args{:} 'psychnoplot' 'cpseq'};
else
    if psych
        args = {args{:} 'psych'};
    end
    if seq == 1
        args = {args{:} 'seqt'};
    elseif seq == 2
        args = {args{:} 'sequence'};
    end
end
if DATA.plot.showcp == 1
    args = {args{:} 'cpt'};
elseif DATA.plot.showcp == 2
    args = {args{:} 'cp'};
elseif DATA.plot.showcp == 3
    args = {args{:} 'cphist'};
end

if DATA.plot.showem
    args = {args{:} 'eyem'};
end
if DATA.state.uselfp
    args = {args{:} 'lfpt'};
end

if DATA.plot.nmin > 0
    args = {args{:} 'nmin' DATA.plot.nmin};
end
if DATA.plot.sdfw > 0 
    args = {args{:} 'sdfw' DATA.plot.sdfw};
end
if strfind(Expt.Header.expname,'tfXip')
    args = {args{:} 'sxcx'};
end

if DATA.plot.collapse
    args = {args{:} 'collapse' 2};
end
if DATA.plot.flip
    args = {args{:} 'flip'};
end
if DATA.plot.acov
    args = {args{:} 'acov'};
end

    
function mousept = myellipse(mousept, finish)


%
%
% mousept modes
%  1 draw ellipse
%  2 re-size both dimensions
%  3 rotate
%  4 move ellipse
%  5 change left edge (move center and h radius)
if nargin > 1
    start = mousept.start;
else
    mousept.mode = 4;
end
if mousept.mode == 1
    a = (finish(1,1)-start(1,1))/2;
    b = (finish(2,2)-start(2,2))/2;
    mousept.r = abs([a b]);
    mousept.c = [finish(1,1)+start(1,1) finish(2,2)+start(2,2)]/2;
elseif mousept.mode == 2
    a = (finish(1,1)-mousept.c(1));
    b = (finish(2,2)-mousept.c(2));
    mousept.r = abs([a b]);
elseif mousept.mode == 3
    
    dy = (finish(2,2)-mousept.c(2))./mousept.yrange;
    dx = (finish(1,1)-mousept.c(1))./mousept.xrange;
    t = -dy/dx;
    mousept.angle = atan(-dy/dx);
    a = mousept.r(1);
    b = mousept.r(2);
%   [a,b] = exy(mousept.angle,mousept.r(1),mousept.r(2));
%    fprintf('%.3f %.3f %.3f\n',mousept.angle,a,b);
%    a = (mousept.r(1) * cos(mousept.angle)) + mousept.ratio * mousept.r(1) * sin(mousept.angle)
%    b = (mousept.r(2) * cos(mousept.angle)) - (mousept.r(2) * sin(mousept.angle)/mousept.ratio)
    %b = mousept.r(2);
elseif mousept.mode == 4  %% move ellipse
    a = mousept.r(1);
    b = mousept.r(2);
elseif mousept.mode == 5  %% move ellipse
    a = mousept.r(1);
    b = mousept.r(2);
    mousept.c(1) = finish(1,1)-mousept.offset(1);
    mousept.c(2) = finish(2,2)-mousept.offset(2);
elseif mousept.mode == 6 %% move R boundary
    a = (finish(1,1)-mousept.c(1));
    mousept.c(1) = mousept.c(1) + (a-mousept.r(1))/2;
    mousept.r(1) = abs(finish(1,1) - mousept.c(1));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 7 %% move L boundary
    a = abs(mousept.c(1) - finish(1,1));
    mousept.c(1) = mousept.c(1) - (a-mousept.r(1))/2;
    mousept.r(1) = abs(finish(1,1) - mousept.c(1));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 8 %% move Top boundary
    b = (finish(2,2)-mousept.c(2));
    mousept.c(2) = mousept.c(2) + (b-mousept.r(2))/2;
    mousept.r(2) = abs(finish(2,2) - mousept.c(2));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 9 %% move Bottom boundary
    b = abs(finish(2,2)-mousept.c(2));
    mousept.c(2) = mousept.c(2) - (b-mousept.r(2))/2;
    mousept.r(2) = abs(finish(2,2) - mousept.c(2));
    b = mousept.r(2);
    a = mousept.r(1);

else
    a = mousept.r(1);
    b = mousept.r(2);
end

a = a./mousept.xrange;
b = b./mousept.yrange;
sn = sin(mousept.angle);
cn = cos(mousept.angle);
%sn = sin(pi/4);
%cn = cos(pi/4);
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = mousept.xrange .* (x .* cn + y .*sn) + mousept.c(1);
yr = mousept.yrange .* (y .* cn - x .*sn) + mousept.c(2);
mousept.lasth = plot(real(xr),real(yr),'color',mousept.color);

function KeyPressed(a, ks)
global mousept;

get(a);
mousept.mode
if strmatch(ks.Key,'delete') & mousept.mode == 5
    DeleteCluster(mousept.cluster, a);
elseif ks.Key == 'n'
    NewCluster(a);
end

function NewCluster(a)
DATA = GetDataFromFig(a);
j = length(DATA.cluster);
while isempty(DATA.cluster{j}) & j > 1
    j = j-1;
end
newc = j+1;
it = findobj('Tag','Clusterid');
set(it,'value',newc);

function DeleteCluster(cl,callfig)
DATA = GetDataFromFig(callfig);
if ~isempty(DATA.cluster{cl}) & isfield(DATA.cluster{cl},'h') &  ishandle(DATA.cluster{cl}.h)
    delete(DATA.cluster{cl}.h);
end
DATA.cluster{cl} = [];
if isfield(DATA.Expts{DATA.currentexpt}.gui,'spkrange')
    ispk = DATA.Expts{DATA.currentexpt}.gui.spkrange;
    ispk = [ispk(1):ispk(2)];
    spks = find(DATA.AllData.Spikes.codes(ispk,2) == cl);
    DATA.AllData.Spikes.codes(ispk(spks),2) = 0;
    PlotSpikeXY(DATA,ispk(spks),DATA.spkcolor{1});
% Shouldn't need this. not deleting all clusters. And if we were would
% imply that we are setting not clusters
% DATA.Expts{DATA.currentexpt}.gui.classified = 0;
end
set(DATA.toplevel,'UserData',DATA);
if cl > 1
    newc = cl-1;
elseif length(DATA.cluster) > cl & ~isempty(DATA.cluster{cl+1})
    newc = cl+1
else
    newc = 1;
end
it = findobj('Tag','Clusterid');
set(it,'value',newc);

function FigButtonReleased(src, data)
global mousept;
mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
mousept.down = 0;
mousept.lasth = 0;
DATA = GetDataFromFig(src);
ExcludeTrials(DATA, mousept);

function ExcludeTrials(DATA, mousept)
    if DATA.state.plotseq == 2
        id = find([DATA.Expt.Trials.Trial] > mousept.start(1,1) & ...
            [DATA.Expt.Trials.Trial] < mousept.finish(1,1));
        ex = [DATA.Expt.Trials(id).Trial];
        if isfield(DATA.Expt,'ExcludeCluster')
            DATA.Expt.ExcludeCluster{1} = union(DATA.Expt.ExcludeCluster{1}, ex);
        else
            DATA.Expt.ExcludeCluster{1} = ex;
        end
    end

        
    hold off;    
    DATA.Expt = PlotCombined(DATA, DATA.Expt);
    set(DATA.toplevel,'UserData',DATA);

function mousept = myrect(mousept, finish)
    start = mousept.start;
    mousept.finish = finish;
    mousept.siz = finish-start;
    if mousept.mode == 1 % horizontal line
        mousept.lasth = plot([start(1,1) finish(1,1)],[start(2,2) start(2,2)],'r');
    else
        mousept.lasth = plot([start(1,1) finish(1,1)],[start(2,2) finish(2,2)],'r');
    end

function FigButtonPressed(src, data)
global mousept;
set(src, 'WindowButtonMotionFcn',@FigButtonDragged);

mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
mousept.down = 1;
mousept.lasth = 0;
mousept.drags = 0;
hold on; %othewise drawing ellipse deletes data
if mousept.mode == 1
    mousept.start = get(gca,'CurrentPoint');
end
mousept.l = mousept.start;
mousept.siz = [0 0];

function FigButtonDragged(src, data)
global mousept;

if mousept.down
    pt = get(gca,'CurrentPoint');
    if mousept.lasth
        delete(mousept.lasth);
    end
    mousept= myrect(mousept,pt);
%   plot(pt(1,1),pt(2,2),'+');
    mousept.drags = mousept.drags+1;
end

function ButtonPressed(src, data)
global mousept;

DATA = GetDataFromFig(src);

it = findobj('Tag','Clusterid');
if ~isempty(it)
    mousept.cluster = get(it,'value');
else
    mousept.cluster = 1;
end

mousept.drags = 0;

mousept.color = DATA.spkcolor{mousept.cluster+1};
mousept.lasth = 0;
lastkey = get(gcf,'CurrentCharacter');

% don't want the axis rescaling as we draw the ellipse
set(gca,'Xlimmode','Manual','Ylimmode','Manual');

oldmode = mousept.mode;
mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
mousept.down = 1;
if mousept.cluster <= length(DATA.cluster)
    C = DATA.cluster{mousept.cluster};
else
    C.h = 0;
end
hold on; %othewise drawing ellipse deletes data
if mousept.mode == 1
    mousept.start = get(gca,'CurrentPoint');

    if length(DATA.cluster) >= mousept.cluster & ...
            isfield(DATA.cluster{mousept.cluster},'h') & ...
            ishandle(DATA.cluster{mousept.cluster}.h)
        for j = 1:length(DATA.cluster)
            distance(j) = DistanceToCluster(DATA.cluster{j},mousept.start(1,1:2));
        end
        if(min(distance) > 1.05) % pressed outside = start over
            delete(C.h);
        else  %% pressed inside; select this cluster, move if mouse moves
            [d, cl]= min(distance);
            mousept.cluster = cl;
            C = DATA.cluster{cl};
            if ~isempty(it)
                set(it,'value',cl);
            end
            DATA.cluster{cl}.h = DrawCluster(DATA.cluster{cl}, DATA.spkcolor{cl+1});
            mousept.mode = 5;
            mousept.c = [DATA.cluster{cl}.x(1) DATA.cluster{cl}.y(1)];
            mousept.r = [DATA.cluster{cl}.x(2) DATA.cluster{cl}.y(2)];
            mousept.offset = mousept.start(1,1:2) - mousept.c;
            set(DATA.cluster{cl}.h,'linewidth',2);
            if oldmode == 5 %second press in ellipse - move it
                mousept.down = 1;
                mousept.lasth = DATA.cluster{cl}.h;
            else
                mousept.down = 0;  %%ignore drag, release
            end
            set(DATA.toplevel,'UserData',DATA);
        end
    end
elseif mousept.mode == 2 %R button
     if C.h & ishandle(C.h) 
         delete(C.h); 
     end
     cl = mousept.cluster;
     mousept.c = [DATA.cluster{cl}.x(1) DATA.cluster{cl}.y(1)];
     mousept.r = [DATA.cluster{cl}.x(2) DATA.cluster{cl}.y(2)];
    mousept.start = get(gca,'CurrentPoint');
     if mousept.start(1,1) > mousept.c(1) + mousept.r(1)
         mousept.mode = 6;
     elseif mousept.start(1,1) < mousept.c(1) - mousept.r(1)
         mousept.mode = 7;
     elseif mousept.start(2,2) > mousept.c(2) + mousept.r(2)
         mousept.mode = 8;
     elseif mousept.start(2,2) < mousept.c(2) - mousept.r(2)
         mousept.mode = 9;
     end
     fprintf('Start %.2f,%.2f C %.2f,%.2f, R%.2f %.2f, mode %d\n',...
         mousept.start(1,1),mousept.start(2,2),mousept.c(1),mousept.c(2),mousept.r(1),mousept.r(2),mousept.mode);
elseif mousept.mode == 3 % R button
     if C.h delete(C.h); end

end


function distance = DistanceToCluster(C, pos);
   
if isempty(C) | ~isfield(C,'x');
    distance = NaN;
    return;
end
xy = pos - [C.x(1) C.y(1)];
xy = xy ./ [C.x(3) C.y(3)];
cn = cos(-C.angle);
sn = sin(-C.angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = p./[C.x(2)./C.x(3) C.y(2)./C.y(3)];
distance = sum(distance.^2);

function PlotSpikeXY(DATA, spkid, color)
plot(DATA.Spikes.cx(spkid),DATA.Spikes.cy(spkid),'.',...
    'color',color,'markersize',DATA.ptsize);

function ClassifySpikes(mousept,src,varargin)

DATA = GetDataFromFig(src);
cl = mousept.cluster;
DATA.currentcluster = cl;
C.x = [mousept.c(1) mousept.r(1) mousept.xrange];
C.y = [mousept.c(2) mousept.r(2) mousept.yrange];
C.angle = -mousept.angle;
C.h = mousept.lasth;
C.params = [DATA.plot.clusterX DATA.plot.clusterY];
C.Arange = DATA.clusterArange;
C.Brange = DATA.clusterBrange;
C.Erange = DATA.clusterErange;
if isfield(DATA,'spklist') && ~isempty(DATA.spklist)
   expspks = DATA.spklist;
else
    expspks = DATA.spkrange(1):DATA.spkrange(2);
end
C.firstspk = expspks(1);
C.lastspk = expspks(end);
if DATA.currenttrial > 1
    DATA.cluster{cl}.Cluster = C;
    DATA.cluster{cl}.lastspk = C.firstspk-1;
else
    DATA.cluster{cl} = C;
end
DATA.newcluster(mousept.cluster) = 1;
%mousept.lasth
colors = mycolors;


DATA = SetSpkCodes(DATA,expspks,1);
DATA.state.recut = 1;
SetGui(DATA);
id = find(DATA.AllData.Spikes.codes(expspks,2) == 0);
if ~DATA.densityplot
    PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{1});
end
if DATA.state.autoplot
    DATA = CountSpikes(DATA, DATA.currentexpt,'replot');
else
    DATA = CountSpikes(DATA, DATA.currentexpt);
end
DATA.Expts{DATA.currentexpt}.gui.classified = 2;
if DATA.plot.showISI
    GetFigure('ISI');
    isis = CalcISI(DATA.Expts{DATA.currentexpt}.Trials);
    id = find(isis < 1000)
    hist(isis(id),100);
end

set(DATA.toplevel,'UserData',DATA);

function dp = CalcClusterdp(DATA, cl)
    expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
    id = find(DATA.AllData.Spikes.codes(expspks,2) == cl-1);
    nid = find(DATA.AllData.Spikes.codes(expspks,2) ~= cl-1);
    rx = std(DATA.Spikes.cx(expspks(id)));
    cx = mean(DATA.Spikes.cx(expspks(id)));
    ry = std(DATA.Spikes.cy(expspks(id)));
    cy = mean(DATA.Spikes.cy(expspks(id)));
    dx = DATA.Spikes.cx(expspks(id))- mean(DATA.Spikes.cx(expspks(id)));
    dy = DATA.Spikes.cy(expspks(id))- mean(DATA.Spikes.cy(expspks(id)));
    dc = abs(dx+i*dy);
    dx = DATA.Spikes.cx(expspks)- mean(DATA.Spikes.cx(expspks(id)));
    dy = DATA.Spikes.cy(expspks)- mean(DATA.Spikes.cy(expspks(id)));
    d = abs(dx+i*dy);
    [s, did] = sort(d);
    dc = d;
    dc(id) = 0;
    d(nid) = 0;
% need to order these somehow first before doing fpos/drate
% ? rank order d before setting id, nid to zero.

%detection rate
    drate = cumsum(dc(did)) ./ sum(dc);
%false positive rate
fpos = cumsum(d(did)) ./ sum(d);
%area gives ROC
dp = trapz(fpos,drate);


function DATA = DrawXYPlot(DATA, expspks)
    
    ho = ishold;
    
    if DATA.state.recut
        ctype = 2;
    else
        ctype = 1;
    end

    expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
    for j = 1:DATA.nclusters+1
        id = find(DATA.AllData.Spikes.codes(expspks,ctype) == j-1);
        PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{j});
        hold on;
    end
    hold on;
    DATA = DrawClusters(DATA, 0);
        CalcClusterdp(DATA,1);
    FinishXYPlot(DATA);
    if ~ho
        hold off;
    end
    
function FinishXYPlot(DATA)
    if DATA.plot.autoscale == 0
        set(gca,'Xlim',DATA.plot.clusterXrange,'Ylim',DATA.plot.clusterYrange);
    end
    xlabel(DATA.spkvarnames{DATA.plot.clusterX});
    ylabel(DATA.spkvarnames{DATA.plot.clusterY});
    it = findobj(DATA.xyfig, 'Tag','Clusterid');
    c = get(it,'value');
    if ~isempty(DATA.explabels{DATA.currentexpt})
        expname = DATA.explabels{DATA.currentexpt};
%     else
        expname = DATA.Expts{DATA.currentexpt}.Header.expname;
    end
    if DATA.currenttrial > 1
        expname = [expname sprintf('from %d',DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial)];
    end
    if isfield(DATA,'cluster') & ~isempty(DATA.cluster) & length(DATA.cluster) >= c & ~isempty(DATA.cluster{c}) & isfield(DATA.cluster{c},'params')
    title(sprintf('%s: C%d %sX%s',expname,c,DATA.spkvarnames{DATA.cluster{c}.params(1)},DATA.spkvarnames{DATA.cluster{c}.params(2)})); 
    else
    title(sprintf('%s: C%d',expname,c)); 
    end
        

function DATA = SetSpkCodes(DATA, expspks, show)

DATA.AllData.Spikes.codes(expspks,2) = 0;
%
%really want to limit this to spike in the current scope;id = find(DATA.Spks.cluster == 0);

if ~isfield(DATA,'cluster')
 nclusters = 0;
 DATA.cluster = {};
end
nclusters = length(DATA.cluster);
cspks = expspks;
for cl = nclusters:-1:1
    C = DATA.cluster{cl};
    while ~isempty(C) & isfield(C,'x')
        if ~isfield(C,'lastspk')  & ~isempty(DATA.spklist)
            C.lastspk = DATA.spklist(end);
        end
        if ~isfield(C,'firstspk') & ~isempty(DATA.spklist)
            C.firstspk = DATA.spklist(1);
        elseif isfield(C,'firstspk') & C.firstspk > 0
            cspks = expspks(find(expspks >= C.firstspk & expspks <= C.lastspk));
        end
    x = (DATA.Spikes.cx(cspks) - C.x(1))./C.x(3);
    y = (DATA.Spikes.cy(cspks) - C.y(1))./C.y(3);
    xr = x .* cos(C.angle) + y .* sin(C.angle);
    yr = y .* cos(C.angle) - x .* sin(C.angle);
    d = (yr./C.y(2)*C.y(3)).^2 + (xr./C.x(2)*C.x(3)).^2;
    id = cspks(find(d < 1));
    DATA.AllData.Spikes.codes(id,2) = cl;
    if isfield(C,'Cluster')
        C = C.Cluster;
    else
        C = {};
    end
    if ~DATA.densityplot & show
        plot(DATA.Spikes.cx(id),DATA.Spikes.cy(id),'.',...
            'color',DATA.spkcolor{cl+1},'markersize',DATA.ptsize);
        hold on;
    end
    end
end


function ButtonReleased(src, data)
global mousept;

%mousept
if mousept.down == 0 %do nothing
   % delete(mousept.lasth);
   % mousept.lasth = 0;
    return;
end
mousept.down = 0;
pt = get(gca,'CurrentPoint');
    if mousept.lasth & ishandle(mousept.lasth)
        delete(mousept.lasth);
    end
  mousept= myellipse(mousept,pt);
  if mousept.mode == 5
      set(mousept.lasth,'linewidth',2);
  end
ClassifySpikes(mousept, src);
%mousept.drags

function ButtonDragged(src, data)

global mousept;

if mousept.down
    pt = get(gca,'CurrentPoint');
    if mousept.lasth && ishandle(mousept.lasth)
        delete(mousept.lasth);
    end
    mousept= myellipse(mousept,pt);
%    plot(pt(1,1),pt(2,2),'+');
    mousept.drags = mousept.drags+1;
end


function Expts = CheckSaccades(Expts, emfile)

load(emfile);
saccth = 10;

for j = 1:length(Expts)
    crate = Expt.Header.CRrates(1) * 10000;
    for k = 1:length(Expts{j}.Trials)
        if strfind(Expts{j}.Trials(k).OptionCode,'+2a')
            bs = Expts{j}.Trials(k).TrialStart;
            es = Expts{j}.Trials(k).End(end);
            [tdiff, emt] = min(abs([Expt.Trials.Start] - bs));
            if tdiff < 160
                maxi = min([length(Expt.Trials(emt).lh) length(Expt.Trials(emt).rh) length(Expt.Trials(emt).lv) length(Expt.Trials(emt).rv)]);
                ch = (Expt.Trials(emt).lh(1:maxi)+Expt.Trials(emt).rh(1:maxi))/2;
                cv = (Expt.Trials(emt).lv(1:maxi)+Expt.Trials(emt).rv(1:maxi))/2;
                dv = smooth(diff(cv),5);
                dh = smooth(diff(ch),5);
                v = 10000 * sqrt(dv.^2+dh.^2)./crate;
                first = (bs -Expt.Trials(emt).ftime)/crate;
                last =  (es -Expt.Trials(emt).ftime)/crate;
                times = ftime + [1:length(cv)] .*crate;
                saccs = find(v > 10);
                id = find(saccs > last);
                if ~isempty(id)
                    [rate, pk] = max(v(saccs(id)));
                    pk = id(pk);
                    sv(1) = mean(cv(first:saccs(pk)-10));
                    sv(2) = mean(cv(saccs(pk)+10:end));
                    sh(1) = mean(ch(first:saccs(pk)-10));
                    sh(2) = mean(ch(saccs(pk)+10:end));
                    saccsize = abs(diff(sh) + i* diff(sv));
                    sacdir = atan2(diff(sv),diff(sh));
                end
            end
        end
    end
end


function name = BuildName(name)

if isempty([strfind(name,'/') strfind(name,'\')])
    name = [pwd '/' name];
end
id = strfind(name,':');
if id
    if isunix
        name = ['/bgc' name(id(1)+1:end)];
        name = strrep(name,'\','/');
    else
        name = name(id(1)+1:end);
    end
else
    name = name;
end
