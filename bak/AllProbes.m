function result = AllProbes(varargin)
%
% Make sure lfpblank resp is included when doing all files.
% RLS OT probe 1 needs scaling.
% Plot for LFP timecouse averaged over all stimuli, by layer
% add y, n to exps
name = 'AllProbes';
%
%
%If strings is empty, no list is shown. Otherwise a listbox is included.
%
strings = {};
AllExpts = {};
DATA = [];
result = [];
na = 1;
toplevel = 0;
if  isstruct(varargin{1}) & isfield(varargin{1},'Title')
    tag = [name varargin{1}.Title];
    name = tag;
elseif  isstruct(varargin{1}) & isfield(varargin{1},'toplevel')
    DATA = varargin{1};
    toplevel = DATA.toplevel;
    na = 2;
elseif isstruct(varargin{1})
    AllExpts = varargin{1};
elseif length(varargin) >2 & strcmp(varargin{2},'Tag')
    tag = varargin{3};
elseif length(varargin) >2 & ishandle(varargin{1}) & isfigure(varargin{3})
    toplevel = varargin{3};
    tag = get(toplevel,'Tag');
    na = 4;
elseif length(varargin) >1 & isfigure(varargin{1}) 
    toplevel = varargin{1};
    tag = get(toplevel,'Tag');
    na = 2;
elseif ischar(varargin{1}) && strcmp(varargin{1},'list')
    tag = varargin{2};
elseif ischar(varargin{1}) && exist(varargin{1},'file')
tag = varargin{1}; 
else
tag = name; 
end
init = 0;
if toplevel == 0 
toplevel = findobj('Tag',tag);
end

if (isempty(toplevel) || toplevel == 0)%% called by user for new window
    if ischar(varargin{1}) & strcmp(varargin{1},'list') % a list of expt files
        AllExpts = ReadExptList(varargin{2});
        nprobes = AllExpts.nprobes;
    elseif ischar(name) & isempty(strfind(varargin{1},'.mat')) & exist(varargin{1},'file')
        fid = fopen(varargin{1},'r');
        names = textscan(fid,'%s');
        fclose(fid);
        k = 1;
        for j = 1:length(names{1})
            path = names{1}{j};
            if exist(path,'file')
                load(path);
                nprobes(k) = size(AllExpts.exps{1}.lfp,4);
                AllExpts.file = path;
                AllCells{k} = AllExpts;
                AllCells{k}.Track = LoadDriftFile(strrep(path,'.all.mat','.drift.mat'));
                k = k + 1;
            else 
                fprintf('No File %s\n',path);
            end
        end            
    elseif iscell(varargin{1}) % a list of AllExpt files
        names = varargin{1};
        k = 1;
        for j = 2:length(names)
            [pdir, monk] = fileparts(names{1});
            if isempty(monk)
                if strfind(names{1},'lem')
                    monk = 'lem';
                end
            end
            path = [names{1} '/' names{j} '/' monk names{j} '.all.mat'];
            if exist(path,'file')
                load(path);
                nprobes(k) = size(AllExpts.exps{1}.lfp,4);
                AllExpts.file = path;
                AllExpts.type = 0;
                AllCells{k} = AllExpts;
                AllCells{k}.Track = LoadDriftFile(strrep(path,'.all.mat','.drift.mat'));
                AllCells{k}.nprobes = nprobes(k);
                k = k + 1;
            else 
                fprintf('No File %s\n',path);
            end
        end
        
    elseif exist(varargin{1},'file')
    [pdir, tag] = fileparts(varargin{1});
    name = tag;
    load(varargin{1});
    AllExpts.file = varargin{1};
    AllExpts.type = 0;
    nprobes(1) = size(AllExpts.exps{1}.lfp,4);
    else
     fprintf('No file %s\n',varargin{1});
     return;
    end
end
if ~isempty(toplevel)
  if strncmpi(varargin{1},'store',5)
    set(topelevel,'UserData',varargin{2});
    DATA = varargin{2};
  elseif isempty(DATA)
    DATA = get(toplevel,'UserData');
  end
end

if ~isempty(AllExpts)
    DATA = SetAllFile(DATA,AllExpts);
end

IsBusy(DATA);

if nargin
    if strncmpi(varargin{na},'update',5)
        update;
    elseif strncmpi(varargin{na},'close',5)
        f = fields(DATA.tag);
        for j = 1:length(f)
            CloseTag(DATA.tag.(f{j}));
        end
    elseif strncmpi(varargin{na},'getstate',5)
        result = DATA;
        
    elseif strncmpi(varargin{na},'PlotLamFind',5) %plots that help identify IV
        PlotExptGroup(varargin{na});
    elseif strncmpi(varargin{na},'plottype',5)
        na = na+1;
        DATA.plot.type = varargin{na};
    elseif strncmpi(varargin{na},'combine',5)
        for k = 1:length(DATA.AllCells)
        nid = 0;
        for j = 1:length(DATA.AllCells{k}.names)
            if strfind(DATA.AllCells{k}.names{j},DATA.plot.compare)
                nid = nid+1;
                ids(nid) = j;
            end
        end
        if nid
            if isfield(DATA.AllCells{k}.Track,'probes')
                probes = DATA.AllCells{k}.Track.probes;
            else
                probes = DATA.plot.plotlfps;
            end
            lfps{k} = DATA.AllCells{k}.exps{ids(1)}.lfp(:,:,:,probes);
            lfpblanks{k} = DATA.AllCells{k}.exps{ids(1)}.lfpblank(:,probes);
            lfptimes{k} = DATA.AllCells{k}.exps{ids(1)}.lfptimes;
            lfpsizes(k,:) = size(lfps{k});
            lfpn{k} = DATA.AllCells{k}.exps{ids(1)}.lfpn;
            lfpx{k} = DATA.AllCells{k}.exps{ids(1)}.x;
            lfpwrs{k} = DATA.AllCells{k}.exps{ids(1)}.lfpwr;
            lfpf{k} = DATA.AllCells{k}.exps{ids(1)}.lfpfrq;
            Headers{k} = DATA.AllCells{k}.exps{ids(1)}.Header;
            Headers{k}.exptype = name2type(DATA.AllCells{k}.names{ids(1)});
        end
        end
        GetFigure(DATA.tag.MainFig);
        subplot(1,1,1);
        type = get(findobj(DATA.toplevel,'Tag','combtype'),'value');
        if type == 3
        shift(1) = 0;
        for j = 2:length(lfps)
        [shift(j), details] = AlignMatrices(lfps{1},lfps{j},[2 3]);
        end
        elseif type == 2
        for j = 1:length(lfps)
            [a,b] = LFPLatency(lfpblanks{j},lfptimes{j});
            [latency(j), shift(j)] = min(a(:,3));
        end
        else
            shift = zeros(1,length(lfps));
        end
        shift = shift - min(shift); 
        [res.lfp, a] = CombineLFPS(lfps,lfpsizes,lfpn, lfpwrs, shift);
        res.lfptimes = mean(cell2mat(lfptimes'));
        res.lfpn = mean(a.nsum,3);
        res.lfpwr = a.plfp;
        res.lfpfrq = lfpf{1};
        res.x = mean(cell2mat(lfpx'));
        res.Header = Headers{1};
        args = PlotArgs(DATA);
        if ndims(squeeze(res.lfp)) ==2
            imagesc(squeeze(res.lfp));
        else
            GetFigure(DATA.tag.MainFig);
            res.Stimvals = DATA.AllCells{k}.exps{ids(1)}.Stimvals;
            PlotAllProbes(res,DATA.plot.type, args{:});
        end
    elseif strncmpi(varargin{na},'compare',5)
        args = PlotArgs(DATA);
        rcs = [];
        for k = 1:length(DATA.AllCells)
        figure(DATA.figures(k));
        nid = 0;
        for j = 1:length(DATA.AllCells{k}.names)
            if strfind(DATA.AllCells{k}.names{j},DATA.plot.compare)
                nid = nid+1;
                ids(nid) = j;
            end
        end
        if nid
            if strmatch(DATA.plot.type,{'AllBlanks' 'FindLam' 'ExpSpikes'})
                DATA.AllCells{k}.probelist = [1:DATA.nprobes];
                DATA.AllCells{k}.exptypes = DATA.AllCells{k}.names;
                PlotExptGroup(DATA.AllCells{k},DATA.plot.type);
            elseif strmatch(DATA.plot.type,{'SpTrigLFP'})
                args = PlotArgs(DATA);
                PlotSpTrigLFP(DATA,ids(1),'cell', k, args{:});
                set(gca,'XTick',[])
            elseif strmatch(DATA.plot.type,{'Psych'})
                PlotPsych(DATA.AllCells{k}.exps{ids(1)})
            else
                rcs{k} = PlotAllProbes(DATA.AllCells{k}.exps{ids(1)},DATA.plot.type,'probes',DATA.probelist, args{:});
                set(gca,'Xtick',[],'Ytick',[]);
                rcs{k}.mapprobe =DATA.AllCells{k}.mapprobe;
            end
        end
        end
        if length(rcs)
            DATA.rcs = rcs;
            savercs = 1;
            if savercs
                save('rcs.mat','rcs');
            end
            GetFigure('Psych');
            PlotRCs(DATA);
        end
        set(DATA.toplevel,'UserData',DATA);    
    elseif strncmpi(varargin{na},'wintile',7)
        ncells = length(DATA.AllCells);
        if ncells > 15
            nr = 4; nc = 5;
        else
            nr = 2; nc = 3;
        end
        mx = get(0,'monitorpositions')
        if size(mx,1) == 1
            xo = 100;
            fw = (mx(3)-xo)./nc;
            fh = (mx(4)-xo)./nr;
        else
            fw = (mx(2,3)-mx(2,1))./nc;
            fh = mx(2,4)./nr;
            xo = mx(2,1);
        end
        ir = 0;
        ic = 0;
        for j = 1:length(DATA.AllCells)
            figpos(j,:) = [(xo +ir * fw ) (ic * fh) fw fh];
            set(DATA.figures(j),'Position',figpos(j,:));
            ic = ic+1;
            if ic >= nc-1
                ic = 0;
                ir = ir + 1;
            end
        end
        scrsz = get(0','screensize');
        
    elseif strncmpi(varargin{na},'recompare',7)
        GetFigure('Psych');
            PlotRCs(DATA);
    elseif strncmpi(varargin{na},'locatecells',7)
        a = GetFigure('AllSummary');
        hold off;
        
        CalcCellDepths(DATA);

    elseif strncmpi(varargin{na},'checklfps',7)
        a = GetFigure('AllSummary');
        hold off;      
        DATA = CheckLFPs(DATA);
    elseif strncmpi(varargin{na},'lfpcohs',6)
        a = GetFigure('AllSummary');
        hold off;      
        DATA = CalcLFPCoherences(DATA);
    elseif strncmpi(varargin{na},'savedrift',5)
        ElectrodeDrift = DATA.Track.drift;
        Track = DATA.Track;
        if isfield(DATA,'Trak')
        Track = DATA.Trak;
        end
        if isfield(DATA,'Markers')
            Track.Markers = DATA.Markers;
        end
        Track.probes = DATA.probelist(DATA.plot.plotlfps);
        driftfile = [DATA.dir '/' strrep(DATA.name,'.all','.drift.mat')];
        save(driftfile,'ElectrodeDrift','Track');
    elseif strncmpi(varargin{na},'loaddrift',5)
        driftfile = [DATA.dir '/' strrep(DATA.name,'.all','.drift.mat')];
        load(driftfile);
        DATA.Track.drift = ElectrodeDrift;
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(varargin{na},'makeplot',5)
     figure(DATA.mfig);
     subplot(1,1,1);
     PlotExptGroup(DATA,DATA.plot.maketype);
    elseif strncmpi(varargin{na},'setentry',5)
     n = get(DATA.listui, 'value');
     figure(DATA.mfig);
     subplot(1,1,1);
     if strmatch('onestim', DATA.plot.type)
         GetFigure(DATA.tag.OneStimFig);
     end
     if strmatch(DATA.plot.type,{'AllBlanks' 'FindLam' 'ExpSpikes'})
         PlotExptGroup(DATA,DATA.plot.type);
     elseif strmatch(DATA.plot.type,{'TrackBlanks'})
         if ~isfield(DATA,'LFP')...
                 || ~strcmp(DATA.LFP.Header.fstring,DATA.fstrings{n})
             name = strrep([DATA.dir '/' DATA.fstrings{n}],'.c1.','.lfp.');
             load(name);
             DATA.LFP = LFP;
             DATA.LFP.Header.fstring = DATA.fstrings{n};
             set(DATA.toplevel,'UserData',DATA);
         end
         if isfield(DATA,'TrackBlank') && ~strcmp(DATA.TrackBlank.name,DATA.fstrings{n})
             rebuild = 1;
         else
             rebuild = 0;
         end
         if ~isfield(DATA,'TrackBlank') || rebuild
             GetFigure('TrackBlanks');
             [a,b] = TrackBlanks(DATA.LFP,DATA.plot.plotlfps);
             DATA.TrackBlank = b;
             DATA.TrackBlank.peaks = a;
             DATA.TrackBlank.name = DATA.fstrings{n};
         else
             a = DATA.TrackBlank.peaks;
             b = DATA.TrackBlank;
         end
         GetFigure('Tracking');
         plot(b.xtimes, a,'o-');
         hold on;
         dv = mean(b.maxs)-mean(b.maxs);
         
         
        if strfind(DATA.fstrings{n},'rds')
            stimname = 'RDS';
            
        elseif strfind(DATA.fstrings{n},'grating')
            stimname = 'Grating';
        elseif strfind(DATA.fstrings{n},'image')
            stimname = 'Image';
        elseif strfind(DATA.fstrings{n},'rls')
            stimname = 'RLS';
        elseif strfind(DATA.fstrings{n},'square')
            stimname = 'Squarewv';
        else
            stimname = 'Other';            
        end
        for j = DATA.exps{n}.Header.expids
            DATA.Track.Stimnames{j} = stimname;
        end

        Track = DATA.Track;
         subplot(2,1,1);
         hold off;
         Track.BlankMax(DATA.exps{n}.Header.expids) = b.maxs;
         [l,p] = min(b.latencies');
         Track.BlankLat(DATA.exps{n}.Header.expids) = p;
         Track.xtimes(DATA.exps{n}.Header.expids) = b.xtimes;
         varmax = b.maxvs;
         [maxv, maxvp] = max(varmax');
         Track.VarMax(DATA.exps{n}.Header.expids) = maxvp;
         Track.frameresp(DATA.exps{n}.Header.expids,:) = b.frameresp;
         [l,p] = min(b.vlatencies');
         Track.VarLat(DATA.exps{n}.Header.expids) = p;
         id = find(b.maxvart) < 1000;
         [maxr, maxt] = max(b.maxvs);
         DATA.Track = Track;
         PlotDrift(DATA,b);
         set(DATA.toplevel,'UserData',DATA);
     elseif strmatch(DATA.plot.type,{'SpTrigLFP'})
         args = PlotArgs(DATA);
         PlotSpTrigLFP(DATA,n,args{:});
     elseif strmatch(DATA.plot.type,{'CellTrigLFP'})
         args = PlotArgs(DATA);
         PlotSpTrigLFP(DATA,n,args{:},'cells');
     else
         hold off;
         args = PlotArgs(DATA);
         DATA.exps{n}.Header.exptype = DATA.exptypes{n};
         if strcmp(DATA.plot.type,'Test')
         PlotAllProbes(DATA.exps{n},'LFPeigvectors', args{:});
         elseif strcmp(DATA.plot.type,'Psych')
             PlotPsych(DATA.exps{n},args{:});
         elseif strcmp(DATA.plot.type,'CollapseX')
         rc = PlotAllProbes(DATA.exps{n},'onestim', args{:},'sumx');
         else
         rc = PlotAllProbes(DATA.exps{n},DATA.plot.type, args{:});
         end
     end
    elseif strncmpi(varargin{na},'store',5)
        set(DATA.toplevel,'UserData',DATA);
    else
        j = 1;
        init = 1;
        while(j < nargin)
            if(strncmpi(varargin{j},'name',3))
                j = j+1;
                name = varargin{j};
            end
            j = j+1;
        end
    end
else
    init = 1;
end

if init & isempty(findobj('Tag',tag))
  DATA.name = name;
  DATA.tag.top = tag;
  if exist('AllCells','var')
      for j = 1:length(AllCells)
          AllCells{j}.mapprobe = zeros(1,nprobes(j));
          if strfind(AllCells{j}.file,'lemM040')
              AllCells{j}.mapprobe(8) = 7;
          end
      end
      DATA.AllCells = AllCells;
  end
  DATA.nprobes = max(nprobes);
  InitInterface(DATA, DATA.fstrings, DATA.suffs);
else
    NotBusy(DATA);
end


function AllExpts = ReadExptList(name)

fid = fopen(name,'r');
a = textscan(fid,'%s');
names = a{1};

for j = 1:length(names);
    AllExpts.names{j} = names{j};
    AllExpts.exps{j} = PlotAllProbes(names{j});
    nprobes(j) =   size(AllExpts.exps{j}.lfp,4);
    AllExpt.type = 0;
end
AllExpts.file = name;
AllExpts.type = 1;
AllExpts.nprobes = nprobes;

function PlotPsych(rc, varargin)
        
                id = find(rc.pk.xv > -1000);
        if isfield(rc,'pp')
            subplot(2,1,2);
            hold off;
            nmin = max([rc.pp.n])./5;
            fitpsf(rc.pp,'nmin',nmin, 'showfit','shown');
            subplot(2,1,1);
        else
            subplot(1,1,1);
        end
        hold off;
        plot(rc.pk.xv(id),rc.pk.kernel(id));
        rc.pk.xv = rc.pk.xv(id);
        rc.pk.kernel = rc.pk.kernel(id);
        rc.pk.n = rc.pk.n(id,:);
        sigs = unique(rc.or);
        sigs= mod(sigs,180);
        for j = 1:length(sigs)
            hold on;
            plot([sigs(j) sigs(j)],[0 max(rc.pk.kernel)],'r:');
        end
        title(sprintf('mean %.1f frames',mean(rc.pk.n(:))));

function PlotRCs(DATA,varargin)


    GetFigure(DATA.tag.MainFig);
    if strcmp(DATA.plot.compare,'DCOR')
        CompareDCOR(DATA,varargin);
        GetFigure('Latencies');
%        R = CompareLatencies(DATA,'plot');
    elseif strcmp(DATA.plot.compare,'OPRC') || strcmp(DATA.plot.compare,'PPRC')
        GetFigure('Latencies');
        hold off;
        R = CompareLatencies(DATA,'plot');
    end
    
function CompareOPRC(DATA,varargin)
    R = CompareLatencies(DATA,'plot');
    
function R = CompareLatencies(DATA,varargin)
    plottype = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'plot',3)
            plottype = 1;
        end
        j = j+1;
    end
    
    rcs = DATA.rcs;
    k = 0;
    for j = 1:length(rcs)
        if isfield(rcs{j},'latencies')
            k = k+1;
        R.slat(k,:) = rcs{j}.latencies;
        R.ssig(k,:) = rcs{j}.latsig;
        R.llat(k,:) = rcs{j}.lfplatencies(:,1)';
        R.lsig(k,:) = rcs{j}.lfplatsig;
        R.cells(k,:) = rcs{j}.iscell;
        end
    end
    for j = 1:size(R.slat,2)
        
       id = find(R.ssig(:,j) > 10 & R.lsig(:,j) > 10 & ~isnan(R.llat(:,j)));
       R.mumean(j) = mean(R.slat(id,j));
       R.mun(j) = length(id);
       R.lmean(j) = mean(R.llat(id,j));
       cid = find(R.ssig(:,j) > 10 & R.lsig(:,j) > 10 & ~isnan(R.llat(:,j)) & R.cells(:,j) == 1);
       R.sun(j) = length(cid);
       R.sumean(j) = mean(R.slat(cid,j));
    end
    
    if plottype == 1
    plot(R.sumean);
    hold on; 
    plot(R.lmean,'r');
    plot(R.mumean,'g');
    end
        
function CompareDCOR(DATA,varargin)
    rcs = DATA.rcs;
    cellcp = []; cellpcp = []; celllcp = [];
    cellpos = [];
    cellcv = [];
    cellid = [];
    cellacp = [];
    abscps = [];
    probes = [];
    cps = [];
    cellpwrcp = [];
    alln = zeros(100,1);
    allbands = zeros(4,100,4);
    for j = 1:length(rcs)
        pk(j,:) = rcs{j}.pk.kernel;
        pkx(j,:) = rcs{j}.pk.xv;
        upstim(j) = mean([rcs{j}.cp.dnstim]);
        probes(j,:) = (rcs{j}.probes - rcs{j}.blankmax).* DATA.Track.probesep;
        cp = [rcs{j}.cp.cp];
        id = find(sum(diff(cat(1,rcs{j}.cp.sigcount),[],2),2) < 0);
        cp(id) = 1-cp(id);
        %cps is MU CP by counts on signal trials
        cps(j,:) = cp;
        %abscps ignores stim prefs
        abscps(j,:) = [rcs{j}.cp.cp];
        %pcps uses pref dir of subspace map
        pcps(j,:) = [rcs{j}.cp.pcp]; %from subspace pref
        cvs(j,:) = rcs{j}.cv;
        allcp(j,:) = [rcs{j}.cp.cp];

        if isfield(rcs{j},'bandsigpwr')
            if DATA.plot.cptype == 5
                pwrsign = sum(rcs{j}.bandsigpwr(:,3:4),2);
            elseif DATA.plot.cptype == 4
                pwrsign = sum(rcs{j}.bandsigpwr(:,1:2),2);
            else
                pwrsign = sum(rcs{j}.bandsigpwr,2);
            end
            cp = [rcs{j}.cp.cp];
            id = find(pwrsign > 0);
            cp(id) = 1-cp(id);
            lpwrcps(j,:) = cp;
        end
        if isfield(rcs{j},'lfpprefdir')
            cellproj(:,1) = cos(rcs{j}.lfpprefdir-mean([rcs{j}.cp.upstim] * pi/180));
            cellproj(:,2) = cos(rcs{j}.lfpprefdir-mean([rcs{j}.cp.dnstim] * pi/180));
            cp = [rcs{j}.cp.cp];
            id = find(diff(cellproj,[],2) > 0);
            cp(id) = 1-cp(id);
            lcps(j,:) = cp;
        end

        
        if isfield(rcs{j},'cellcp')
            cid = [];
            for k = 1:length(rcs{j}.spkres)
                cid(k) = ~isempty(rcs{j}.spkres{k});
            end
            cid = find(cid);
            for ic = 1:length(cid)
            %cp not adjusted for pref
            cellacp = [cellacp rcs{j}.spkres{cid(ic)}.cp.cp];
            end
            %cellcp pref is assigned from counts from signal trials.
            cellcp = [cellcp rcs{j}.cellcp(cid)];
            %cp based on pref of the LFP response in teh subspace map
            celllcp = [celllcp rcs{j}.celllcp(cid)];

            %cp based on pref of the LFP response in teh subspace map
            if isfield(rcs{j},'bandsigpwr')
                for ic = 1:length(cid)
                    spk = rcs{j}.spkres{cid(ic)};
                    probe = round(mean(spk.probes));
                    if rcs{j}.mapprobe(probe)
                        probe = rcs{j}.mapprobe(probe);
                    end
                    if pwrsign(probe) > 0
                    cellpwrcp = [cellpwrcp 1-spk.cp.cp];
                    else
                    cellpwrcp = [cellpwrcp spk.cp.cp];
                    end
                end
                        
            end
            %cp based on the cells own subspace map. From direction of mean
            %netspk vector.
            cellpcp = [cellpcp rcs{j}.cellpcp(cid)];
            cellpos = [cellpos (rcs{j}.cellprobe(cid) - rcs{j}.blankmax).*DATA.Track.probesep];
            cellcv = [cellcv rcs{j}.cellcv(cid)];
            cellid = [cellid ones(size(rcs{j}.cellcv(cid))) .*j];
            rid(length(cid)) = 0;
            rid(find(rid == 0)) = j;
        end
        if strcmp(DATA.plot.type,'LFPpwrdiff')
            iprobes(j,:) = 20 +rcs{j}.probes - rcs{j}.blankmax;
            allbands(:, iprobes(j,:),1,j) = rcs{j}.pwrbands;
            s =  sign(rcs{j}.bandsigpwr');
            allbands(:, iprobes(j,:),2,j) = rcs{j}.pwrbands .*s;
            s = abs(sin(rcs{j}.x(rcs{j}.prefs,1) .*pi/180));
            s = s > sin(pi/4);
            s = repmat(sign(s' -0.5),4,1);
            allbands(:, iprobes(j,:),3,j) = rcs{j}.pwrbands .*s;
            s = abs(sin(rcs{j}.prefdir .*pi/180));
            s = s > sin(pi/4);
            s = repmat(sign(s -0.5),4,1);
            allbands(:, iprobes(j,:),4,j) = rcs{j}.pwrbands .*s;
            alln(iprobes(j,:),:) = alln(iprobes(j,:),:) + 1 ;
        end
    end
    subplot(3,1,1);
    hold off;
    [a,b] = sort(upstim);
    imagesc([min(pkx(:)) max(pkx(:))], 1:length(b), pk(b,:));
    hold on;
    plot(upstim(b),1:length(b),'w:');
    colorbar;
    subplot(3,1,3);
    hold off;
    if ~isfield(DATA.plot,'muCVcrit')
        DATA.plot.muCVcrit = 0.2;
    end
    if ~isfield(DATA.plot,'cellCVcrit')
        DATA.plot.cellCVcrit = 0.7;
    end
    muid = find(cvs > DATA.plot.muCVcrit);
    id = find(cellcv >= DATA.plot.cellCVcrit);
    uid = find(cellcv < DATA.plot.cellCVcrit);
    if DATA.plot.cptype == 1 %by pref from subspace map
        mucps = pcps;
        cellcps = cellpcp;
    elseif DATA.plot.cptype == 0 %by counts on signal trials
        mucps = cps;
        cellcps = cellcp;
    elseif DATA.plot.cptype ==2 %by lfp subspace pref
        mucps = lcps;
        cellcps = celllcp;
    elseif ismember(DATA.plot.cptype, [3 4 5]) %LFP pref by power in bands for sig trials
        mucps = lpwrcps;
        cellcps = cellpwrcp;
    else %abs CP for upstim(first), ignoring all prefs
        mucps = abscps;
        cellcps = cellacp;
    end

    plot(probes(muid), mucps(muid),'o');
    hold on;
    plot(cellpos(id),cellcps(id),'ko','markerfacecolor','k');
    plot(cellpos(uid),cellcps(uid),'ro','markerfacecolor','r');
    title(sprintf('Mean CP %.3f, Cells %.3f(%.3f)n=%d(%d)',mean(mucps(muid)),...
            mean(cellcps(id)),mean(cellcps),length(id),length(cellcps)));

    plot([min(probes(muid)) max(probes(muid))],[0.5 0.5],'k:');
    % allbands(,probe,cptype,cell)
    % cptype = 1 : power in this band
    % cptype = 2 : power in this band
    if strcmp(DATA.plot.type,'LFPpwrdiff')
        GetFigure('LFP power by choice');
        id = find(alln > 0);
        nm = repmat(alln(id)',[4 1 4]);
        allsd = squeeze(std(allbands(:,id,:,:),[],4))./sqrt(nm);
        mband = squeeze(sum(allbands(:,id,:,:),4))./nm;
        for j = 1:4
            subplot(4,1,j);
            title('Pref from power in band');
            errorbar(squeeze(mband(:,:,j))',squeeze(allsd(:,:,j))');
        end
    end



function PlotDrift(DATA, b)
    
    Track = DATA.Track;
    for j = 1:length(Track.Stimnames)
        if isempty(Track.Stimnames{j})
            stimtypes(j) = 0;
        else
            stimtypes(j) = strmatch(Track.Stimnames{j},{'RDS', 'Grating', 'Image', 'RLS', 'Squarewv', 'Other'});
        end
    end
    id = find(Track.BlankMax > 0 & stimtypes > 1)  %all expts where this has been done
    if nargin > 1
        hold off;
    plot(b.xtimes,b.maxs,'mo-'); %Max blank resp
    hold on;
    else
        hold off;
    end
    if length(id)
        c = (mean(Track.frameresp(id,2))+mean(Track.BlankMax(id)))/2;
        plot(Track.xtimes(id),Track.drift(id)+c,'g');
        hold on;
        plot(Track.xtimes(id),Track.BlankMax(id),'r');
        plot(Track.xtimes(id),Track.VarMax(id),'b');
        plot(Track.xtimes(id),Track.frameresp(id,1),'k');
        plot(Track.xtimes(id),Track.frameresp(id,2),'k:');
        plot(Track.xtimes(id),Track.frameresp(id,3),'k:');
        plot(Track.xtimes(id),Track.BlankLat(id),'r:');
        plot(Track.xtimes(id),Track.VarLat(id),'b:');
        legend('drift','blankmax','varmax','Frame')
    else
        c = 0;
    end
    if nargin > 1
    plot(b.xtimes,b.blankshifts+c,'ro-'); %from xcorr
    hold on;
    plot(b.xtimes,b.varshifts+c,'bo-'); %from xcorr
    subplot(2,1,2);
    hold off;
    imagesc(b.tslice);
    hold on;
    plot(b.maxs,1:size(b.tslice,1),'w:');
    end
  
function type = name2type(name)
        if strfind(name, 'OTRC');
            type = 'OTRC';
        elseif strfind(name, 'PPRC');
            type = 'PPRC';
        elseif strfind(name, 'OPRC');
            type = 'OPRC';
        elseif strfind(name, 'ODRC');
            type = 'ODRC';
        elseif strfind(name, 'square.CO');
            type = 'Flash';
        elseif strfind(name, 'DCORRC');
            type = 'DCORRC';
        else
            id = regexp(name,'[A-Z][A-Z][A-Z]*');
            if isempty(id)
            type = 'Unknown';
            else
            type = name(id(1):end);
            end
        end

        
function IsBusy(DATA)
if isfield(DATA,'toplevel')
    set(DATA.toplevel,'Name','Busy...');
end

function NotBusy(DATA)
if isfigure(DATA.toplevel)
    set(DATA.toplevel,'Name',DATA.tag.top);
end

function DATA = SetAllFile(DATA,AllExpts)
    DATA.exps = AllExpts.exps;
    DATA.fstrings = AllExpts.names;
    [pdir, name] = fileparts(AllExpts.file);

    DATA.dir = pdir;
    for j = 1:length(DATA.exps)
        DATA.suffs{j} = '';
        if isfield(DATA.exps{j},'Header')
            ts = DATA.exps{j}.Header.BlockStart(1)./(60 * 10000); % time in min
            DATA.suffs{j} = sprintf('ed%.2f,%.1f',mean(DATA.exps{j}.Header.eds),ts);
        end
        DATA.exptypes{j} = name2type(DATA.fstrings{j});
        DATA.exps{j}.Header.exptype = DATA.exptypes{j};
    end
    if AllExpts.type == 1
        return;
    end
    driftfile = [DATA.dir '/' strrep(name,'.all','.drift.mat')];
    DATA.Track = LoadDriftFile(driftfile);
    if isfield(DATA.Track,'probes')
        DATA.plot.plotlfps = DATA.Track.probes;
    end
    DATA = SetExptDrifts(DATA);
    
    function Track = LoadDriftFile(driftfile)
        
Track = [];
if exist(driftfile,'file')
        load(driftfile);
        Track.drift = ElectrodeDrift;
        if exist('Track','var')
            Track = Track;
            if ~isfield(Track,'drift')
                Track.drift = ElectrodeDrift;
            end
        end
    end
    if ~isfield(Track,'probesep')
        Track.probesep = 50;
    end



    function DATA = SetExptDrifts(DATA)
% add drift values to Headers of each expt
nexp = length(DATA.exps);
nb = 0;
oldn = 1;
for j = 1:nexp
    starts = DATA.exps{j}.Header.BlockStart;
    nb = nb+length(starts);
    allstarts(oldn:nb) = starts;
    expid(oldn:nb) = j;
    blkid(oldn:nb) = 1:length(starts);
    oldn = nb+1;
end
[a,b] = sort(allstarts);
for j = 1:length(b);
    ei = expid(b(j));
    H = DATA.exps{ei}.Header;
    if isfield(DATA,'Track') && isfield(DATA.Track,'drift')
    DATA.exps{ei}.Header.drift(blkid(b(j))) = DATA.Track.drift(j);
    end
    DATA.exps{ei}.Header.expids(blkid(b(j))) = j;
    if isfield(H,'BlockTrial')
    DATA.Track.BlockStart(j) = DATA.exps{ei}.Header.BlockTrial(blkid(b(j)));
    end
end


function [mlfp, details] = CombineLFPS(lfps, sz, lfpn, lfpwrs,shift)

nr = lfpn;
for j = 1:length(lfpwrs)
    nfreqs(j) = size(lfpwrs{j},1);
end
nfreq = mean(nfreqs);
if sum(std(sz)) == 0 %same size
    ntimes = mean(sz(:,1));
   nsum = zeros(size(nr{1}));
   np = max(sz(:,4));
   nsum = repmat(nsum,[1 1 np+max(shift)]);
    for j = 1:length(lfps)

        pwr = lfpwrs{j};
        for p = 1:min(sz(:,4)) % smallest # of probes
            slfp (:,:,:,p+shift(j),j) = lfps{j}(:,:,:,p) .* repmat(reshape(nr{j},[1 size(nr{j})]),[ntimes 1 1]);
            if ndims(lfpwrs{j}) == 3
                slfpwr(:,:,p+shift(j),j) = pwr(:,:,p) .* sum(nr{j}(:));
            else
                slfpwr(:,p+shift(j),j) = pwr(:,p) .* sum(nr{j}(:));
            end
            nsum(:,:,p+shift(j)) = nsum(:,:,p+shift(j))+nr{j};
        end
    end
    details.nsum = nsum;
    psum = squeeze(sum(nsum,2));
%        mlfp(:,:,:,p) = sum(slfp(:,:,:,p,:),5) ./ repmat(reshape(nsum,[1 size(nsum)]),[ntimes 1 1]);
   for p = 1:size(nsum,3)
       mlfp(:,:,:,p) = sum(slfp(:,:,:,p,:),5) ./ shiftdim(repmat(nsum(:,:,p),[1 1 ntimes]),2);
            if ndims(lfpwrs{j}) == 3
                details.plfp(:,:,:,p) = sum(slfpwr(:,:,p,:),4) ./ sum(psum(:,p));
            else
                details.plfp(:,:,:,p) = sum(slfpwr(:,p,:),3) ./ sum(psum(:,p));
            end
   end
%    slfp = sum(slfp,5) ./ repmat(reshape(nsum,[1 size(nsum)]),[ntimes 1 1]);
%        slfp = sum(cat(5,lfps{:}),5);
else
    fprintf('Dimensions don''t match:\n');
    for j = 1:size(sz,1)
        fprintf('%d ',sz(j,:));
        fprintf('\n');
    end
    mlfp = [];
    return;
end


function PlotSpTrigLFP(DATA, n, varargin)

[lp, it] =  GetPop('plotstyle',DATA.toplevel);
str = get(it,'String')
DATA.plot.style = str(lp,:);
lp = lp-1;
diffy =  GetCheck('DiffY',DATA.toplevel);
sumy =  GetCheck('SumY',DATA.toplevel);
meanp =  GetCheck('MeanProbe',DATA.toplevel);
celln = 0;
Expt = [];
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'cell',4)
        if isfield(DATA,'AllCells')
        j = j+1;
        celln = varargin{j};
        Expt = DATA.AllCells{celln}.exps{n};
        else
            celln = 1:length(DATA.exps{n}.spkres);
        end
    end
    j = j+1;
end
if isempty(Expt)
    Expt = DATA.exps{n};
end

if ~isfield(Expt,'sptriglfp')
    dfile = [];
    if isfield(Expt.Header,'filename')
        dfile = Expt.Header.filename;
    end
    if isempty(dfile) | ~exist(dfile,'file')
        dfile = [DATA.dir '/' DATA.fstrings{n}];
    end
    sfile = regexprep(dfile,'.c1.*mat','.sptrig.mat');
    load(sfile);
    if celln & isfield(DATA,'AllCells')
    for j = 1:length(DATA.AllCells{celln}.exps)
        DATA.AllCells{celln}.exps{j}.sptriglfp = sptrig{j}.lfp;
        DATA.AllCells{celln}.exps{j}.sptriglfpn = sptrig{j}.lfpn;
        DATA.AllCells{celln}.exps{j}.triglfptimes = sptrig{j}.lfptimes;
        if isfield(sptrig{j},'cells')
            DATA.AllCells{celln}.exps{j}.trigcells = sptrig{j}.cells;
        end
    end
        Expt = DATA.AllCells{celln}.exps{n};
    else
    for j = 1:length(DATA.exps)
        DATA.exps{j}.sptriglfp = sptrig{j}.lfp;
        DATA.exps{j}.sptriglfpn = sptrig{j}.lfpn;
        DATA.exps{j}.triglfptimes = sptrig{j}.lfptimes;
        if isfield(sptrig{j},'cells')
            DATA.exps{j}.trigcells = sptrig{j}.cells;
        end
    end
    Expt = DATA.exps{n};
    end
end
spv = Expt.sptriglfp;
iscell = zeros(1,size(spv,2));
if isfield(Expt,'trigcells')
    if celln
        args = {'smooth', [0 DATA.plot.smoothk(2)]};
        if isfield(Expt.Header,'LFPamps')
           % args = {args{:} 'scales' Expt.Header.LFPamps};
        end
        if strcmp(DATA.plot.style,'Pcolor')
            args = {args{:} 'avg'};
        end
        if DATA.plot.diffy
            args = {args{:} 'diffy'};
        end
        trig{1}.name  = Expt.Header.filename;
        trig{1}.cells = Expt.trigcells;
        TrigCSD(trig,'cell',celln,'autooff',args{:});
        return;
    end
    cellprobes = [Expt.trigcells.probe];
    for j = 1:length(cellprobes)
        if size(Expt.trigcells(j).lfptrig,1) == size(spv,1)
        spv(:,:,:,round(cellprobes(j))) = Expt.trigcells(j).lfptrig;
        iscell(round(cellprobes(j))) = 1;
        end
    end
else
    cellprobes = [];
end
np = size(spv,2);
[nr,nc] = Nsubplots(np);
cr = [min(spv(:)) max(spv(:))];
plotlfps = DATA.plot.plotlfps;
for j =1:np
    if iscell(j) % not just a cell, but has been put in above
        utype = 'S';
    else
        utype = 'M';
    end
    subplot(nr,nc,j);
    if sumy
        
    [a, delays(j,:)] = min(squeeze(sum(spv(50:150,plotlfps,:,j),3)));
    end
    if diffy & size(Expt.sptriglfpn,2) > 1
        nsp = sprintf('%.0f-%.0f',Expt.sptriglfpn(j,1),Expt.sptriglfpn(j,2));
    elseif sumy & size(Expt.sptriglfpn,2) > 1
        nsp = sprintf('%d+%d',Expt.sptriglfpn(j,1),Expt.sptriglfpn(j,2));
    else 
        nsp = sprintf('%.0f',Expt.sptriglfpn(j,1));
    end
    tiv = Expt.triglfptimes;
    if lp
        if size(spv,3) > 1 && diffy
            plot(tiv,squeeze(spv(:,plotlfps,1,j) - spv(:,plotlfps,2,j)));
        elseif size(spv,3) > 1 && sumy
            plot(tiv,squeeze(spv(:,plotlfps,1,j) + spv(:,plotlfps,2,j))/2);
        else
            if meanp
                plot(tiv,squeeze(mean(spv(:,plotlfps,1,j),4)));
            else
                plot(tiv,squeeze(spv(:,plotlfps,1,j)));
            end
            if size(spv,3) > 1
                hold on;
                if meanp
                plot(tiv,squeeze(mean(spv(:,plotlfps,2,j),2)),'r');
                else
                    plot(tiv,squeeze(spv(:,plotlfps,2,j)),':');
               end
                
                hold off;
            end
                
        end
    else
        if size(spv,3) > 1 && diffy
            hold off;
            imagesc(squeeze(spv(:,:,1,j) - spv(:,:,2,j))');
        elseif size(spv,3) > 1 && sumy
            hold off;
            imagesc(tiv, 1:np,squeeze(spv(:,:,1,j) + spv(:,:,2,j)/2)');
            hold on; plot(delays(j,:)-50,plotlfps,'w')
        else
            if DATA.plot.yvals(1) <= size(spv,3) & DATA.plot.yvals(1) > 0
            imagesc(squeeze(spv(:,:,DATA.plot.yvals(1),j))');
            else
            imagesc(squeeze(spv(:,:,1,j))');
            end
        end
    end
    h = title(sprintf('P%d %s %cU',j,nsp,utype));
    if utype == 'S';
        set(h,'color','r');
    end
end
if lp == 0
    colorbar;
    subplot(nr,nc,j+1);
    if sumy
    imagesc(delays);
    end
end



function PlotExptGroup(DATA, type)

RDS = 1;
ofig = gcf;
probes = DATA.plot.plotlfps;
if strncmpi(type,'FindLam',7)
    id = strmatch('Flash',DATA.exptypes);

    nr = 1;
    if ~isempty(id)
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'default','timerange',[400 1000 500 1500]);
        results{nr}.title = 'Flash';
        nr = nr+1;
    end
    id = strmatch('PPRC',DATA.exptypes);
    if ~isempty(id)
        res = PlotAllProbes(DATA.exps{id(1)},DATA.plot.type,'LFPeig');
        pp=sum(sum(res.eigresp,3),2);
        if ~isempty(res.eigblank) && sum(res.eigblank) < 0
            [a, xi] = max(pp);
        else
            [a, xi] = min(pp);
        end
        [a, ui] = max(pp);
        [a, li] = min(pp);
% find extremum that is nearest the middle of the range  
        if abs(ui - length(pp)/2) > abs(li - length(pp)/2)
            xi = li;
        else
            xi = ui;
        end
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'onestim','xvals',xi,'timerange',[400 1000 500 1500]);
%        latency = LFPLatency(results{nr}.lfpim.z,results{nr}.lfptimes);
        results{nr}.title = sprintf('PP = %.2f',res.x(xi));
        nr = nr+1;
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'VarBlank','timerange',[400 1000 500 1500]);
        results{nr}.title = sprintf('PPVAR',res.x(xi));
        nr = nr+1;
    end
    id = strmatch('OPRC',DATA.exptypes);
    if ~isempty(id)
        res = PlotAllProbes(DATA.exps{id(1)},DATA.plot.type,'LFPeig');
        pp=sum(sum(res.eigresp,3),2);
        if ~isempty(res.eigblank) && sum(res.eigblank) < 0
            [a, xi] = max(pp);
        else
            [a, xi] = min(pp);
        end
        [a, ui] = max(pp);
        [a, li] = min(pp);
% find extremum that is nearest the middle of the range  
        if abs(ui - length(pp)/2) > abs(li - length(pp)/2)
            xi = li;
        else
            xi = ui;
        end
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'onestim','xvals',xi,'timerange',[400 1000 500 1500]);
        results{nr}.title = sprintf('OP = %.2f',res.x(xi));
        nr = nr+1;
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'VarBlank','timerange',[400 1000 500 1500]);
        results{nr}.title = sprintf('OPVAR = %.2f',res.x(xi));
        nr = nr+1;
    end
    id = strmatch('DCORRC',DATA.exptypes);
    if ~isempty(id)
        res = PlotAllProbes(DATA.exps{id(1)},DATA.plot.type,'LFPeig');
        pp=sum(sum(res.eigresp,3),2);
        if ~isempty(res.eigblank) && sum(res.eigblank) < 0
            [a, xi] = max(pp);
        else
            [a, xi] = min(pp);
        end
        [a, ui] = max(pp);
        [a, li] = min(pp);
% find extremum that is nearest the middle of the range  
        if abs(ui - length(pp)/2) > abs(li - length(pp)/2)
            xi = li;
        else
            xi = ui;
        end
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'onestim','xvals',xi,'timerange',[400 1000 500 1500]);
        results{nr}.title = sprintf('DCOR = %.2f',res.x(xi));
        nr = nr+1;
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'VarBlank','timerange',[400 1000 500 1500]);
        results{nr}.title = sprintf('DCORV');
        nr = nr+1;
    end
    id = strmatch('OTRC',DATA.exptypes);
    if ~isempty(id)
        results{nr} = PlotAllProbes(DATA.exps{id(1)},'blank','timerange',[400 1000 500 1500]);
        results{nr}.title = sprintf('OTRC Blank');
        nr = nr+1;
    end
    nr = nr-1;
    for j = 1:nr
        rc = results{j};
        mp = NaN;
        lp = NaN;
        subplot(1,nr,j);
        hold off;
        imagesc(results{j}.lfpim.x,results{j}.lfpim.y,results{j}.lfpim.z);
        hold on;
        if isfield(results{j},'lfplatencies')
            if ndims(results{j}.lfplatencies) == 3
            latency = squeeze(results{j}.lfplatencies(1,1,:))./10;
            else
            latency = results{j}.lfplatencies(:,3)./10;
            end
            plot(latency,results{j}.lfpim.y,'r:');
            [lat, lp] = min(latency);
        end
        if isfield(rc,'lfpblank')
            
                blanklat = LFPLatency(rc.lfpblank, rc.lfptimes);
        [maxr,maxt] = max(abs(rc.lfpblank));
        [a,blankmaxp] = max(maxr);
        [blankl, blanklp] = min(blanklat(:,3));
        else
            blanklp = 0;
            blankmaxp = 0;
        end

        [maxr, maxt] = max(abs(results{j}.lfpim.z'));
        [a, mp] = max(maxr);
        if isfield(results{j},'lfpzc')
            plot(results{j}.lfpztimes./10,results{j}.lfpzc,'r:');
        end
        if isfield(results{j},'nrep')
            nrep = results{j}.nrep;
        else
            nrep = 0;
        end
        title(sprintf('%s N%.0f\nLat%d(%d) Max%d(%d)',results{j}.title,nrep,lp,blanklp,mp,blankmaxp));
    end
elseif strncmpi(type,'TrackAllBlank',8)
    bid = [];
    for j = 1:length(DATA.exps)
        if isfield(DATA.exps{j},'lfpblank') & ~isempty(DATA.exps{j}.lfpblank)...
                & isfield(DATA.exps{j},'mlfp');
            bid = [bid j];
        end
    end
    DATA.plot.type = 'TrackBlanks';
    set(DATA.toplevel,'UserData',DATA);
    for j = 1:length(bid)
        set(DATA.toplevel,'Name',sprintf('Track:%s',DATA.exptypes{bid(j)}));
        drawnow;
        set(DATA.listui,'value',bid(j));
        AllProbes(DATA.toplevel,'setentry');
        DATA = get(DATA.toplevel,'UserData');
        set(DATA.toplevel,'UserData',DATA);
    end
    
elseif strncmpi(type,'AllBlank',5)
    id = strmatch('Flash',DATA.exptypes);
    nr = 1;
    if ~isempty(id)
         flsh = PlotAllProbes(DATA.exps{id(1)},'default','timerange',[400 1200 500 1500],'csd');
         csdz = mean(flsh.csdzc) - mean(DATA.exps{id(1)}.Header.drift);
    end

    bid = [];
    for j = 1:length(DATA.exps)
        if isfield(DATA.exps{j},'lfpblank') & ~isempty(DATA.exps{j}.lfpblank) & isfield(DATA.exps{j},'mlfp')
            bid = [bid j];
            mins(j) = min(DATA.exps{j}.lfpblank(:));
            maxs(j) = max(DATA.exps{j}.lfpblank(:));
            if strfind(DATA.fstrings{j},'rds')
            stimtypes(length(bid)) = RDS;
            else
            stimtypes(length(bid)) = 0;
            end
        end
    end
    cr = [min(mins(bid)) max(maxs(bid))];
    nb = length(bid);
    [r,c] = Nsubplots(nb);
%+v drift means that same signal is on larger probe.
%so subtract drift to keep real position contsant
    for j = 1:length(bid)
        set(0,'currentfig',ofig);
        subplot(1,nb,j);
        drift(j) = mean(DATA.exps{bid(j)}.Header.drift);
        cres = PlotAllProbes(DATA.exps{bid(j)},'onestim','xvals',0,'yvals',0,'csd','timerange',[500 1000],'probes',probes);
        csdmax(j) = cres.csdmax-drift(j);
        csdzc(j,:) = cres.csdzc - drift(j);
        res = PlotAllProbes(DATA.exps{bid(j)},'blank','timerange',[500 1000],'probes',probes);
        vres = PlotAllProbes(DATA.exps{bid(j)},'VarBlank','timerange',[400 1000 500 1500],'subplots',[2 nb j j+nb],'probes',probes);
        vres.blankvar = max(var(vres.lfpim.zb'));
        blankratio(j,1) = vres.blankvar./max(vres.lfpim.z(:));
        blankratio(j,2) = max(vres.blanksdr);
        edepth = mean(res.Header.eds);
        start = res.Header.BlockStart(1)./(60 * 10000);
        if isfield(res,'nrep')
            nrep = res.nrep;
        else
            nrep = 0;
        end
        [maxr,maxt] = max(res.lfpblank);
        [a,maxp] = max(maxr);
        [minl,minp] = min(res.lfplatencies(:,3));
        title(sprintf('%s N%.0f ed%.2f %.1fmin\nLat%.1fP%.0f Max%.0f',DATA.exptypes{bid(j)}, nrep,edepth,start,minl./10,minp,maxp));
        latprobes(j,1) = vres.minlatprobe-drift(j);
        latprobes(j,2) = res.blank.minlat(2)-drift(j);
        maxprobes(j,1) = vres.maxvarprobe-drift(j);
        maxprobes(j,2) = res.blank.maxprobe-drift(j);
        reversals(j,1) = vres.framereverseprobe-drift(j);
        reversals(j,2) = vres.frameminprobe-drift(j);
        reversals(j,3) = vres.framemaxprobe-drift(j);
        nstim(j) = mean(res.lfpn(:));
        caxis(cr);
    end
    GetFigure('AllSummary');
    subplot(1,1,1);
    sid = find(stimtypes ~= RDS & nstim > DATA.plot.rcnmin);
    hold off;
    plot([1:length(sid)],maxprobes(sid,1),'o'); %var max
    hold on;
    plot([1:length(sid)]+0.05,maxprobes(sid,2),'ro'); %blank max
    plot([1:length(sid)],latprobes(sid,1),'s'); %var latency
    plot([1:length(sid)]+0.05,latprobes(sid,2),'rs'); %blank latecny
    plot([1:length(sid)]-0.05,reversals(sid,2),'gs'); %min frameresp
    plot([1:length(sid)],reversals(sid,1),'go','markerfacecolor','g'); %frame reversal
    plot([1:length(sid)]-0.05,reversals(sid,3),'gv'); %max frame resp
    plot([1:length(sid)]+0.1,csdmax(sid),'mo'); 
    plot([1:length(sid)]+0.1,csdzc(sid,1),'m^'); 
    plot([1:length(sid)]+0.1,csdzc(sid,2),'mv'); 
    legend('VarMax','BlnkMax','VarLat','BlnkLat','FrameMin','Reverse','Max','CSD');
    lm = WeightedSum(latprobes(sid,:),nstim(sid));
    mm = WeightedSum(maxprobes(sid,:),nstim(sid));
    fm = WeightedSum(reversals(sid,:),nstim(sid));
    cm = WeightedSum(csdzc(sid,:),nstim(sid));
    cmax = WeightedSum(csdmax(sid),nstim(sid));
    plot([1 length(sid)],[lm(1) lm(1)],'b');
    text(1.1,lm(1)+0.1,'VarLat','color','b');
    plot([1 length(sid)],[lm(2) lm(2)],'r');
    text(1.1,lm(2)+0.1,'BlankLat','color','r');
    plot([1 length(sid)],[mm(1) mm(1)],'b');
    text(1.1,mm(1)+0.1,'MaxVar','color','b');
    plot([1 length(sid)],[mm(2) mm(2)],'r');
    text(1.1,mm(2)+0.1,'BlankMax','color','r');
    plot([1 length(sid)],[csdz csdz],'k');
    set(gca,'xlim',[0.7 length(sid)+0.2]);
    ymin = get(gca,'ylim');
    for j = 1:length(sid)
        text(j-0.1,ymin(1),DATA.fstrings{bid(sid(j))},'rotation',90);
    end
    DATA.Markers.Latency = lm;
    DATA.Markers.Max = mm;
    DATA.Markers.Frame = fm;
    DATA.Markers.csd = [cmax cm];
    DATA.Markers.flashcsdz = csdz;
    set(DATA.toplevel,'UserData',DATA);
 
elseif strncmpi(type,'ExpSpikes',5)
   DATA =  PlotAllExptSpikes(DATA);
   set(DATA.toplevel,'UserData',DATA);
%    PlotExptSpikes(DATA);
elseif strncmpi(type,'Drift',4)
    PlotDrift(DATA);
end

function PlotExptSpikes(DATA)

nexp = length(DATA.exps);
[nr,nc] = Nsubplots(nexp);dbdown
for j = 1:nexp
    subplot(nr,nc,j);
    C = [];
    for probe = 1:length(DATA.probelist);
        ne = length(DATA.exps{j}.Cluster(probe).dprime);
        C(1:ne,probe) = DATA.exps{j}.Cluster(probe).dprime;
        auto(probe) = mean(DATA.exps{j}.Cluster(probe).autocut);
    end
    C(end+1,:) = auto;
    imagesc(C);
end

function DATA = PlotAllExptSpikes(DATA)

ofig = gcf;
subplot(2,1,1);
nexp = length(DATA.exps);
nb = 0;
oldn = 1;
for j = 1:nexp
    starts = DATA.exps{j}.Header.BlockStart;
    nb = nb+length(starts);
    allstarts(oldn:nb) = starts;
    expid(oldn:nb) = j;
    blkid(oldn:nb) = 1:length(starts);
    oldn = nb+1;
end
[a,b] = sort(allstarts);
    C = [];
    H = [];
probes = DATA.plot.plotlfps;
%probes = 10:18;  
for j = 1:length(b);
    ei = expid(b(j));
    Cluster = DATA.exps{expid(b(j))}.Cluster;
    for probe = probes
        if blkid(b(j)) > length(Cluster(probe).dprime)
            C(j,probe) = NaN;
            auto(j,probe) = NaN;
        else
        C(j,probe) = Cluster(probe).dprime(blkid(b(j)));
        auto(j,probe) = DATA.exps{expid(b(j))}.Cluster(probe).autocut(blkid(b(j)));
        if isfield(Cluster,'SpkStats') & ~isempty(Cluster(probe).SpkStats)
        H(j,probe) = Cluster(probe).SpkStats.h(blkid(b(j)),3);
        else
            H(j,probe) = 0;
        end
%        H(j,probe) = Cluster(probe).SpkStats.k(blkid(b(j)));
        end
    end
end


[a, maxe] = max(max(H'));
for j = 1:size(H,1) %# of expts
    cc(j,:) = xcorr(H(j,probes),H(maxe,probes),10,'unbiased');
end
cc = cc ./max(cc(:));
GetFigure('CCAlign');
imagesc(cc);
figure(ofig);
[maxc, drift] =  max(cc');
if ~isfield(DATA.Track,'drift')
DATA.Track.drift = drift-min(drift);
end
for j = 1:length(b);
    ei = expid(b(j));
    DATA.exps{ei}.Header.drift(blkid(b(j))) = DATA.Track.drift(j);
end
subplot(2,2,1);
imagesc(C);
subplot(2,2,3);
hold off;
imagesc(H);
hold on;
plot(drift-10,1:size(H,1),'r:');
plot(DATA.Track.drift+1,1:size(H,1),'w:');
subplot(2,2,[2 4]);
imagesc(auto);
for j = 1:length(b);
    ei  = expid(b(j));
    if isfield(DATA.exps{ei}.Header,'depths')
        text(1,j,[DATA.exptypes{ei} sprintf(' %.2f',DATA.exps{ei}.Header.depths(blkid(b(j))))...
            sprintf(' %.1f',DATA.exps{ei}.Header.BlockStart(blkid(b(j)))/600000)]);
    else
        text(1,j,[DATA.exptypes{expid(b(j))} sprintf(' %.2f',DATA.exps{ei}.Header.eds)]);
    end
end

function args = PlotArgs(DATA)

p = GetPop('OneProbe',DATA.toplevel);
if p > 1
    args = {'probes' p-1};
else
    args = {'probes' DATA.plot.plotlfps};
end
args = {args{:} 'timerange' [DATA.plot.timerange(1) DATA.plot.timerange(2) DATA.plot.timerange(3) DATA.plot.timerange(4)]};
args = {args{:} 'pause' DATA.plot.timerange(5)};
args = {args{:} 'yvals' DATA.plot.yvals};
args = {args{:} 'xvals' DATA.plot.xvals};
args = {args{:} 'freqs' DATA.plot.freqrange};
args = {args{:} 'setsign' DATA.plot.setsign};

if strmatch('onestim', DATA.plot.type)
    args = {args{:} 'figs' {DATA.tag.OneStimFig DATA.tag.MainFig}};
end
if GetCheck('AddBlank',DATA.toplevel);
    args = {args{:} '+blank'};
end
if GetCheck('Interp',DATA.toplevel);
    args = {args{:} 'interpolate'};
end
if GetCheck('ScatterPlot',DATA.toplevel);
    args = {args{:} 'scatterplot'};
elseif GetPop('plotstyle',DATA.toplevel) ==2
    args = {args{:} 'lineplot'};
elseif GetPop('plotstyle',DATA.toplevel) ==3
    args = {args{:} 'csdplot' DATA.plot.smoothk};
end
if GetCheck('SumY',DATA.toplevel);
    args = {args{:} 'sumy'};
end
if GetCheck('DiffY',DATA.toplevel);
    args = {args{:} 'diffy'};
end
if GetCheck('dVdt',DATA.toplevel);
    args = {args{:} 'dvdt'};
end

function SetProbe(a,b)

DATA = GetDataFromFig(a);
probe = get(a, 'value');
if probe == 1
    DATA.probelist = 1:DATA.nprobes;
elseif probe == 26
    DATA.probelist = 1:23;
else
    DATA.probelist = probe-1;
end
set(DATA.toplevel,'UserData',DATA);
AllProbes(DATA.toplevel,'setentry');

function SetPlot(a,b)
DATA = GetDataFromFig(a);
id = get(a, 'value');
strs = get(a,'string');
if strmatch(get(a,'Tag'),'combtype')
    DATA.plot.combtype = strrep(deblank(strs(id,:)),' ','');
elseif strmatch(get(a,'Tag'),'setsign')
    DATA.plot.setsign = id-1;
else
DATA.plot.type = strrep(deblank(strs(id,:)),' ','');
end
set(DATA.toplevel,'UserData',DATA);
AllProbes(DATA.toplevel,'setentry');

    
function MakePlot(a,b)
DATA = GetDataFromFig(a);
id = get(a, 'value');
strs = get(a,'string');
DATA.plot.maketype = strrep(deblank(strs(id,:)),' ','');
set(DATA.toplevel,'UserData',DATA);
AllProbes(DATA.toplevel,'makeplot',DATA.plot.maketype);

function SetComparison(a,b)
DATA = GetDataFromFig(a);
id = get(a, 'value');
strs = get(a,'string');
DATA.plot.compare = strrep(deblank(strs(id,:)),' ','');
set(DATA.toplevel,'UserData',DATA);
AllProbes(DATA.toplevel,'compare');


function Update(a,b,varargin)
DATA = GetDataFromFig(a);
DATA.plot.timerange(1) = GetField('FirstTime',DATA.toplevel);
DATA.plot.timerange(2) = GetField('LastTime',DATA.toplevel);
DATA.plot.timerange(5) = GetField('DelayTime',DATA.toplevel);
DATA.plot.yvals = GetField('Yvals',DATA.toplevel);
DATA.plot.xvals = GetField('Xvals',DATA.toplevel);
DATA.plot.plotlfps = GetField('ShowLFPs',DATA.toplevel);
DATA.plot.autoplot = GetCheck('AutoPlot',DATA.toplevel);
DATA.plot.sumy = GetCheck('SumY',DATA.toplevel);
DATA.plot.diffy = GetCheck('DiffY',DATA.toplevel);
DATA.plot.smoothk = GetField('SmoothKernel',DATA.toplevel);
DATA.plot.rcnmin = GetField('RCNmin',DATA.toplevel);
DATA.plot.freqrange(1) = GetField('minFreq',DATA.toplevel);
DATA.plot.freqrange(2) = GetField('maxFreq',DATA.toplevel);
set(DATA.toplevel,'UserData',DATA);
if DATA.plot.autoplot
    AllProbes(DATA.toplevel,'setentry');
end

function CloseTag(tag)
it = findobj('Tag',tag);
if ~isempty(it)
    close(it);
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

function [value, it] = GetPop(tag, varargin)

if nargin == 2 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    value = get(it(1),'value');
    it = it(1);
else
    value = 0;
    it = 0;
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


function DATA = GetDataFromFig(a)

    DATA = get(a,'UserData');

    if isempty(DATA)
        b = get(a,'parent');
        DATA = get(b,'UserData');
    end
    if isfield(DATA,'parentfigtag')
        DATA = get(findobj('Tag',DATA.parentfigtag),'UserData');
    end
    
function InitInterface(DATA, strings, suffs)

    scrsz = get(0,'Screensize');
    wsc = scrsz(3) /1280;  %scale everything as if 1280 wide
    if scrsz(3) > 1600
        wsc = 1.5;
    end
    size(1) = 380 * wsc;
    size(2) = 300 * wsc;
    if ~isfield(DATA,'nprobes')
        DATA.nprobes = 24;
    end
    DATA.tag.list = [DATA.tag.top 'List'];
   DATA.parentfigtag = DATA.tag.top;
   DATA.probelist = 1:DATA.nprobes;
   DATA.plot.timerange = [20 2000 500 1500 0.01];
   DATA.plot.smoothk = [0.1 2]; %[ time  probe]
   if ~isfield(DATA.plot,'plotlfps')
   DATA.plot.plotlfps = [1:DATA.nprobes];
   end
   DATA.plot.yvals = 1;
   DATA.plot.xvals = 0;
   DATA.plot.type = 'default';
    DATA.plot.sumy = 1;
    DATA.plot.diffy = 0;
    DATA.plot.meanprobe = 0;
    DATA.plot.dvdt = 0;
    DATA.plot.compare = 'DCOR';
    DATA.plot.rcnmin = 500;
    DATA.plot.freqrange = [1 100];
    DATA.plot.interpolate = 0;
    DATA.plot.cptype = 0;
    DATA.plot.setsign = 0;
    DATA.plot.cellCVcrit = 0.7;
    DATA.plot.muCVcrit = 0.2;

    
    HSPACE = 3;
    VSPACE = 2;
    name = DATA.name;

    if scrsz(3) > 1600
        cw = 12;
        ch = 21;
        VSPACE = 3;
    else
        cw = 12;
        ch = 20;
    end
    if isfield(DATA,'AllCells')
        nlines = 6;
    else
        nlines = 5;
    end
    cntrl_box = figure('Position', [10 scrsz(4)-(size(2)+40) size(1) size(2)],...
        'NumberTitle', 'off', 'Tag',DATA.tag.top,'Name',DATA.name);
   DATA.toplevel = cntrl_box;
   
   for j = 1:length(strings)
       labels{j} = [strings{j} suffs{j}];
   end
    if( ~isempty(strings))
        DATA.listui = uicontrol(gcf, 'Style','listbox','String', labels,...
            'Callback', {@AllProbes, DATA.toplevel, 'setentry'} ,'Tag',DATA.tag.list,...
            'Position',[10 10 size(1) size(2)-((ch+VSPACE)*nlines)]);
        DATA.fstrings = strings;
        
    end

    
    bp(1) = HSPACE; bp(3) = 25; bp(2) = size(2)-ch; bp(4) = 22;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@AllProbes, DATA.toplevel, 'next'},...
        'String', '>>', 'Position', bp);
    bp(1) = bp(1) + bp(3) + HSPACE;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback',  {@AllProbes, DATA.toplevel, 'prev'},...
        'String', '<<', 'Position', bp);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 15 * cw;
    uicontrol(gcf,'Style', 'edit','String',sprintf('%d ',DATA.plot.plotlfps),'Position', bp,'Tag','ShowLFPs','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);

       bp(1) = bp(1)+bp(3)+HSPACE;
       bp(3) = cw*5;
  uicontrol(gcf,'style','pop','string','CPrate|CPbypref|CPbylfp|CPbyLFPpwr|Bands1,2|Bands3,4|CPabs', ...
           'Callback', @SetPopPlot, 'Tag','popplot',...
        'position',bp,'value',DATA.plot.cptype+1);
       bp(1) = bp(1)+bp(3)+HSPACE;
       bp(3) = cw*5;
  uicontrol(gcf,'style','pop','string','BlankEig|Smoothness|None|Match blankresp', ...
           'Callback', @SetPopPlot, 'Tag','setsign',...
        'position',bp,'value',DATA.plot.setsign+1);
    
 %New row
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 6 * cw;
    uicontrol(gcf,'Style', 'checkbox','String', '+blank', 'Tag', 'AddBlank', 'Callback', @Update,...
        'Position', bp);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 7 * cw;
    uicontrol(gcf,'Style', 'checkbox','String', 'Scatter', 'Tag', 'ScatterPlot', 'Callback', @Update,...
        'Position', bp);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 4 * cw;
    uicontrol(gcf,'Style', 'checkbox','String', 'Auto', 'Tag', 'AutoPlot', 'Callback', @Update,...
        'Position', bp);

    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'SumY', 'Tag', 'SumY', 'Callback', @Update,...
        'Position', bp,'value',DATA.plot.sumy);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'DiffY', 'Tag', 'DiffY', 'Callback', @Update,...
        'Position', bp,'value',DATA.plot.diffy);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'Mean', 'Tag', 'MeanProbe', 'Callback', @Update,...
        'Position', bp,'value',DATA.plot.meanprobe);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'dVdt', 'Tag', 'dVdt', 'Callback', @Update,...
        'Position', bp,'value',DATA.plot.dvdt);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'Interp', 'Tag', 'Interp', 'Callback', @Update,...
        'Position', bp,'value',DATA.plot.interpolate);
    
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 4 * cw;
    uicontrol(gcf,'Style', 'text','String','Plot','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 7*cw;
    uicontrol(gcf,'style','pop','string','Default|LFPcmp|LFPdiff|LFPeigresp|LFPeigvec|RexpXblank|Blank only|LFPpwr|LFPpwrdiff|LPFsigdiff|LFPband|LFPgamma|Lines|FrameResp|netspk|onestim|Collapse X|onestim mu|monocs|XTProbes|XTprobesmu|TrackBlanks|SpTrigLFP|CellTrigLFP|AllBlanks|FindLam|ExpSpikes|StimVar|Latency|VarBlank|LFPTrial|csdsum|Psych|Test', ...
        'Callback', @SetPlot, 'Tag','plottype',...
        'position',bp,'value',1);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 6*cw;
    uicontrol(gcf,'style','pop','string','Pcolor|Lines|CSD|test', ...
        'Callback', @Update, 'Tag','plotstyle',...
        'position',bp,'value',1);
 
    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 6*cw;
    uicontrol(gcf,'style','pop','string','AllBlanks|TrackAllBlanks|ExpSpikes|FindLam|Drift|test', ...
        'Callback', @MakePlot, 'Tag','MakePlot',...
        'position',bp,'value',1);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 5 * cw;
    uicontrol(gcf,'Style', 'text','String','Probe','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 4*cw;
    uicontrol(gcf,'style','pop','string','All|1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|1:23', ...
        'Callback', @SetProbe, 'Tag','OneProbe',...
        'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'edit','String',sprintf('%d',DATA.plot.yvals),'Position', bp,'Tag','Yvals','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);

    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'edit','String',sprintf('%d',DATA.plot.xvals),'Position', bp,'Tag','Xvals','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);
    
    
    
    bp(2) = bp(2) - ch;
    bp(1) = HSPACE;
      bp(3) = cw * 4;
    uicontrol(gcf,'Style', 'text','String','Time','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.1f',DATA.plot.timerange(1)),'Position', bp,'Tag','FirstTime','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);
      bp(3) = cw * 5;
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.1f',DATA.plot.timerange(2)),'Position', bp,'Tag','LastTime','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',DATA.plot.timerange(5)),'Position', bp,'Tag','DelayTime','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);

  bp(1) = bp(1)+bp(3)+HSPACE;
  bp(3) = cw * 4;
  uicontrol(gcf,'Style', 'text','String','Smth','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  bp(3) = cw * 6;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f ',DATA.plot.smoothk),'Position', bp,'Tag','SmoothKernel','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);
    
  bp(1) = bp(1)+bp(3)+HSPACE;
  bp(3) = cw * 4;
  uicontrol(gcf,'Style', 'text','String','Nmin','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.0f ',DATA.plot.rcnmin),'Position', bp,'Tag','RCNmin','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);

    bp(2) = bp(2) - ch;
    bp(1) = HSPACE;
      bp(3) = cw * 4;
  uicontrol(gcf,'Style', 'text','String','Freq','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.1f ',DATA.plot.freqrange(1)),'Position', bp,'Tag','minFreq','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.0f ',DATA.plot.freqrange(2)),'Position', bp,'Tag','maxFreq','Callback', ...
	    @Update,'Backgroundcolor',[1 1 1]);

    if isfield(DATA,'AllCells')
     for j = 1:length(DATA.AllCells)
         [a, name] = fileparts(DATA.AllCells{j}.file);
         DATA.AllCells{j}.name = name;
          DATA.figures(j) = GetFigure([DATA.name name]);
     end
      figure(DATA.toplevel);
       bp(2) = bp(2) - ch;
       bp(1) = HSPACE;
       bp(3) = cw*6;
       uicontrol(gcf,'Style', 'pushbutton', 'Callback',  {@AllProbes, DATA.toplevel, 'Compare'},...
           'String', 'Compare', 'Position', bp);
       bp(1) = bp(1)+bp(3)+HSPACE;
       bp(3) = cw*7;
       uicontrol(gcf,'style','pop','string','DCOR|OPRC|PPRC|OTRC|ACRC|square.CO|AC|SZOB', ...
           'Callback', @SetComparison, 'Tag','cmptype',...
        'position',bp,'value',1);
       bp(1) = bp(1)+bp(3)+HSPACE;
       bp(3) = cw*6;
       uicontrol(gcf,'Style', 'pushbutton', 'Callback',  {@AllProbes, DATA.toplevel, 'Combine'},...
           'String', 'Combine', 'Position', bp);
       bp(1) = bp(1)+bp(3)+HSPACE;
       bp(3) = cw*8;
       uicontrol(gcf,'style','pop','string','NoShift|Vshift(blankmax)|Vshift(xcorr)|Vshift(varmax)PrefShift|CP|Pysch|CPall', ...
           'Callback', @SetPlot, 'Tag','combtype',...
        'position',bp,'value',1);
    for j = 1:length(DATA.AllCells)
        allnames{j} = DATA.AllCells{j}.name;
    end
       bp(1) = bp(1)+bp(3)+HSPACE;
       bp(3) = cw*10;
  uicontrol(gcf,'style','pop','string',allnames, ...
           'Callback', @SetOneCell, 'Tag','onecell',...
        'position',bp,'value',1);

   end
       
    
  hm = uimenu(gcf,'Label','File');
  uimenu(hm,'Label','Close','Callback', {@AllProbes, DATA.toplevel,'close'});
  uimenu(hm,'Label','Save Drift','Callback', {@AllProbes, DATA.toplevel, 'savedrift'});
  uimenu(hm,'Label','Locate Cells','Callback', {@AllProbes, DATA.toplevel, 'locatecells'});
  uimenu(hm,'Label','Check LFPs','Callback', {@AllProbes, DATA.toplevel, 'checklfps'});
  uimenu(hm,'Label','Calc LFP Cohrence','Callback', {@AllProbes, DATA.toplevel, 'lfpcoh'});
  uimenu(hm,'Label','Tile','Callback', {@AllProbes, DATA.toplevel, 'WinTile'});
  uimenu(hm,'Label','Test','Callback', {@AllProbes, DATA.toplevel, 'menutest'});
    set(gcf,'Menubar','none');

  DATA.tag.MainFig = [DATA.name 'Fig'];
  DATA.tag.OneStimFig = [DATA.name 'OneStim'];
  a = GetFigure(DATA.tag.MainFig);
  DATA.mfig = a;

  set(cntrl_box,'UserData',DATA);
  
  function SetPopPlot(a,b)
     DATA = GetDataFromFig(a);
     cid = get(a,'value');
     DATA.plot.cptype = cid-1;
     if isfield(DATA,'rcs')
         GetFigure('Psych');
         PlotRCs(DATA);
     end
     set(DATA.toplevel,'UserData',DATA);
     

  function SetOneCell(a,b)
     DATA = GetDataFromFig(a);
     cid = get(a,'value');
     cell = DATA.AllCells{cid};
     DATA.fstrings = cell.names;
     DATA.name = cell.name;
     DATA = SetAllFile(DATA,DATA.AllCells{cid});
     for j = 1:length(DATA.fstrings)
         labels{j} = [DATA.fstrings{j} DATA.suffs{j}];
     end
     nf = get(DATA.listui,'value');
     if nf > length(labels)
         nf = length(labels);
     end
     set(DATA.listui,'string',labels,'value',nf);
     set(DATA.toplevel,'UserData',DATA);
  
function [x, details] = TrackBlanks(Expt, probes, varargin)
%
% Calculate blank resp for each expt, and meausre
% location of peak resp, location of zero crossing
% beneath this, to track electrode movement.
%
% details.frameresp(1) is inversion of frame resp at top
% details.frameresp(2) is probe with minumum modulation at Frameresp
%     (usually very close to crossover)
% details.frameresp(3) is probe with maximum mod
sid = [];
j = 1;
fig = gcf;
while j <= length(varargin)
    if strnmcpi(varargin{j},'tsamples',4)
        j = j+1;
        sid = varargin{j};
    end
    j = j+1;
end
blocks = Expt.Header.BlockStart;
tn = [Expt.Trials.Trial];
bExpt = Expt;
[nc,nr] = Nsubplots(length(blocks)-1);
nc = length(blocks)-1;
np = 1;
iprobe = [probes(1):0.1:probes(end)];
res = PlotRevCorAny(Expt,'lfp');
lfp = res.sdfs.extras{1}.lfp-res.sdfs.alllfp;
sumblank = lfp';
[a, maxv] = max(lfp,[],2);
[a, maxt] =  max(a);
if isempty(sid)
    sid = [maxt-5 maxt maxt+5];
end
if length(Expt.Header.BlockStart) == 1
    ntrials(1) = length(tn);
    blocktrials{1} = 1:length(Expt.Trials);
else
ntrials = [];
for j = 1:length(Expt.Header.BlockStart)-1
    blocktrials{j} = find(tn >= blocks(j) & tn <= blocks(j+1));
    ntrials(j) = length(blocktrials{j});
end
    blocktrials{j+1} = find(tn >= blocks(end));
ntrials(j+1) = length(blocktrials{j+1});
end
nc = sum(ntrials > 20);
for j = 1:length(ntrials)
    if ntrials(j) > 20
        trials = blocktrials{j};
        bExpt.Trials = Expt.Trials(trials);
        res = PlotRevCorAny(bExpt,'lfp','nmin',10,'noplot');
        set(0,'CurrentFigure',fig); 
        subplot(1,nc,np);
        np = np+1;
        lfp = res.sdfs.extras{1}.lfp'-res.sdfs.alllfp';
        allblank(j,:,:) = lfp;
        for x = size(res.sdfs.lfp,1):-1:1
            for y = size(res.sdfs.lfp,2):-1:1
                if ~isempty(res.sdfs.lfp{x,y})
                lfps(:,x,y,:) = res.sdfs.lfp{x,y};
                end
            end
        end
        fr = GetFrameReversal(res.sdfs.alllfp,res.sdfs.lfptimes, probes);
        [vlatencies, vx] = LFPLatency(lfps, res.sdfs.lfptimes);
        allvar(j,:,:) = vx.var';
        [latencies, b] = LFPLatency(lfp, res.sdfs.lfptimes);
        imagesc(lfp);
        xtimes(j) =Expt.Trials(trials(1)).TrialStart./(10000 * 60);
        title(sprintf('ed %.2f, %.1fmin',Expt.Header.depths(j),xtimes(j)));
        a = res.sdfs.extras{1}.lfp(maxt,:);
        details.tslice(j,:) = a;
        details.latencies(j,:) = latencies(:,3);
        details.vlatencies(j,:) = vlatencies(:,3);
        details.frameresp(j,:) = [fr.revp fr.minprobe fr.maxprobe];
        details.maxvs(j,:) = vx.maxr;
        details.maxvart(j,:) = vx.maxt;
        [a, maxs(j)] =  max(interp1(probes,a,iprobe,'spline'));
        for k = 1:length(sid)
            zc = diff(sign(res.sdfs.extras{1}.lfp(sid(k),:)));
            id = find(zc(1:length(probes)-1) < 0);
            if length(id)
            id = id(end);
            x(j,k) = interp1(res.sdfs.extras{1}.lfp(sid(k),id:id+1),[id id+1],0);
            else
                x(j,k) = 24;
            end
        end
    else
        allvar(j,:,:) = NaN;
        x(j) = NaN;
    end
end


sumvar = squeeze(mean(allvar,1));
if length(probes) == 8
    nch = 8;
else
nch = 23;
end
nblk = size(allblank,1);
for blk = 1:size(allblank,1)
    [shifts(blk), bxcs(:,blk)] = AlignMatrices(squeeze(allblank(blk,1:nch,:)),...
        sumblank(1:nch,:),1);
    [vshifts(blk), vxcs(:,blk)] = AlignMatrices(squeeze(allvar(blk,1:nch,:)),...
        sumvar(1:nch,:),1);
end

details.varshifts = vshifts; %from xc
details.blankshifts = shifts; %from xc
details.blankxc = bxcs; %xc matric
details.varxc = vxcs;
details.xtimes = xtimes;
maxs(find(maxs == 0)) = 1;
 details.maxs = iprobe(maxs); %blank max
 details.allblank = allblank;
 details.allvar = allvar;
 details.sumblank = sumblank;

% details.maxvs = vx.maxr;

function details = GetFrameReverse(mlfp, lfptimes,  probes)

if ndims(mlfp) == 3
    mlfp = squeeze(sum(mlfp,2));
end

for j = probes
    [a,b] = famp(lfptimes, mlfp(:,j),1/166.66);
    phases(j) = angle(b);
    amps(j) = abs(b);
end
pdiff = diff(phases);
[a,b] = max(abs(pdiff));
if a > pi 
    if pdiff(b) < 0
        id = find(phases >= phases(b));
    else
        id = find(phases < phases(b+1));
    end
   phases(id) = phases(id) + 2 * pi;
%    id = find(phases < 0);
end
details.amps = smooth(amps,2,'gauss');
[a,b] = max(details.amps);
details.maxprobe = probes(b);
[a,b] = min(details.amps(1:b));
details.minprobe = probes(b);

[y,x] =smhist(phases);
[a,b] = max(y);
details.phase = x(b); %dominant phase
details.phases = smooth(phases,2,'gauss');
if mean(phases) < details.phase
    sp = details.phase - pi/2; %90 deg ahead - halfway to 180;
    id = find(details.phases > details.phase);
else
    sp = details.phase + pi/2; %90 deg ahead - halfway to 180;
    id = find(details.phases > details.phase);
end
lastid = find(diff(id) > 1); %find break;
if length(lastid)
id = id(1:lastid(1));
end
[a,b] = max(abs(diff(details.phases(id))));
revp = probes(b);

if sp > max(details.phases)
    if length(id) <= 1
        revp = NaN;
    else
    revp = interp1(details.phases(id),probes(id),sp,'pchip','extrap');
    if revp < 0 || revp > max(probes)
        revp = NaN;
    end
    end
else
    revp = interp1(details.phases(id),probes(id),sp);
end
details.revp = revp;

function [probe, drift] = CellDepth(DATA, E, H)
    last = 0;
    for k = 1:length(E.probestep)
        id = find(E.Trials > last & E.Trials < E.probestep(k));
        last = E.probestep(k)-1;
        probes(id) = E.probes(k);
    end
    if isempty(k) || isempty(probes) %no steps
        probe = mean(E.probes); %should all be the same
    else
    probes(last+1:length(E.Trials)) = E.probes(k+1);
    probe = mean(probes(find(probes > 0)));
    end
    id = find(ismember(H.BlockStart,E.Trials));
    id = find(ismember(DATA.Track.BlockStart,H.BlockStart(id)));
    drift = mean(DATA.Track.drift(id));
    if isnan(drift)
        id(1)
    end

function DATA = CalcLFPCoherences(DATA)
    
    if isfield(DATA,'AllCells');
        for k = 1:length(DATA.AllCells);
            nid = 0;
            for j = 1:length(DATA.AllCells{k}.names)
            if strfind(DATA.AllCells{k}.names{j},DATA.plot.compare)
                nid = nid+1;
                ids(nid) = j;
            end
            end

            if nid
                d = fileparts(DATA.AllCells{k}.file);
                lfpfile = [d '/' strrep(DATA.AllCells{k}.names{ids(1)},'.c1.','.lfp.')];
                if exist(lfpfile,'file')
                    load(lfpfile);
                    [C, X] = LFPCohere(LFP);
                    DATA.AllCells{k}.lfpcoh = C;
                    allamps(k,:) = MaxCohere(C,'freqs',1:100);
                else
                    fprintf('Can''t Read %s\n',lfpfile);
                end
            end
        end
    end
        
    
    if isfield(DATA,'AllCells');
        for j = 1:length(DATA.AllCells);
            lfpfile = strrep(DATA.AllCells{j}.file,'.all.mat','.lfp.mat');
        end
        subplot(2,1,1);
       imagesc(allamps);
       subplot(2,1,2);
       plot(std(allamps'));
       DATA.LFPamps = allamps;
    end
    
function DATA = CheckLFPs(DATA)
    if isfield(DATA,'AllCells');
        for j = 1:length(DATA.AllCells);
            lfpfile = strrep(DATA.AllCells{j}.file,'.all.mat','.lfp.mat');
            if exist(lfpfile,'file')
                load(lfpfile);
                amps = LFPGains(LFP);
                DATA.AllCells{j}.lfpgain = amps;
                allamps(j,1:length(amps)) = amps./max(amps);
                stds(j) = std(allamps(j,find(allamps(j,:) > 0.1)));
                means(j) = mean(allamps(j,find(allamps(j,:) > 0.1)));
            else
                fprintf('Can''t Read %s\n',lfpfile);
            end
        end
        subplot(2,1,1);
       imagesc(allamps);
       subplot(2,1,2);
       plot(stds);
       hold on;
       plot(means,'r');
       DATA.LFPamps = allamps;
    end
    
function CalcCellDepths(DATA)
    
    for e = 1:length(DATA.exps)
       ns(e) =  length(DATA.exps{e}.spkres);     
    end
    probes = ones(length(DATA.exps),max(ns)) .* NaN;
    drift = probes;
    for e = 1:length(DATA.exps)
        if isfield(DATA.exps{e},'spkres')
            for j = 1:length(DATA.exps{e}.spkres)
                E = DATA.exps{e}.spkres{j};
                if ~isempty(E)
                k = E.cellnumber;
                H = E.Header;
                n(e,k) = length(E.Trials);
                [probes(e,k) drift(e,k)] = CellDepth(DATA, E,H);
                end
            end
        end
    end
    cpos = probes-drift;
    DATA.Cellposdat = cpos;
    for j = 1:size(cpos,2)
        id = find(~isnan(cpos(:,j)));
    DATA.Cellpos(j) = sum(cpos(id,j).*n(id,j))./sum(n(id,j));
    end
    set(DATA.toplevel,'UserDATA',DATA);
    PlotCellLocations(DATA);
    
    
function PlotCellLocations(DATA)
    plot(DATA.Cellposdat,'o');
    hold on;
    for j = 1:length(DATA.Cellpos)
        labs{j} = sprintf('%d',j);
    end
    legend(labs);
    nc = size(DATA.Cellposdat,1);
  if isfield(DATA,'Markers')
      plot([1 nc],[DATA.Markers.Frame(1),DATA.Markers.Frame(1)],'r'); %min frame
      plot([1 nc],[DATA.Markers.Latency(2),DATA.Markers.Latency(2)],'m'); %blank latency
      plot([1 nc],[DATA.Markers.Latency(1),DATA.Markers.Latency(1)],'m:'); %Var latency
      plot([1 nc],[DATA.Markers.Max(2),DATA.Markers.Max(2)],'g');      
      plot([1 nc],[DATA.Markers.flashcsdz,DATA.Markers.flashcsdz],'k');      
  end