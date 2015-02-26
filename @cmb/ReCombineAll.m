function out = ReCombineAll(a,b, varargin)
%
% go through all Expts. If a combined expt exists, use this list, then
% reccombine

out = {};
j = 1;
lfponly = 0;
listonly = 0;
oneprobe = 0;
onetype = 0;
mkall = 0;
savefile = 1;
docells = 1;
cargs = {};
recombinenames = {};

while j <= length(varargin)
    if iscellstr(varargin{j})
        recombinenames = varargin{j};
    elseif strncmpi(varargin{j},'cells',5) || strncmpi(varargin{j},'mucells',7)
        docells = 1;
        lfponly = -1; %don't make lfp
        if strncmpi(varargin{j},'mucells',7)
            cargs = {cargs{:}, 'mucells'};
        end
    elseif strncmpi(varargin{j},'lfponly',5)
        lfponly = 1;
    elseif strncmpi(varargin{j},'listonly',5)
        listonly = 1;
    elseif strncmpi(varargin{j},'mkall',5)
        mkall = 1;
    elseif strncmpi(varargin{j},'oneprobe',5)
        oneprobe = 1;
        lfponly = -1;
    elseif strncmpi(varargin{j},'ORBW',4)
        onetype = 1;
    elseif strncmpi(varargin{j},'reapply',6)
        cargs = {cargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'spkonly',5)
        lfponly = -1;
    end
    j = j+1;
end



if isstruct(a)
    DATA = a;
else
    DATA = GetDataFromFig(a);
end

if isempty(DATA.Expts)
    acknowledge('No Expts','combine error');
    return;
end

if ~isfield(DATA,'Templates')
    load('StdTemplate.mat');
    DATA.Templates = Templates;
    DATA.TemplateInfo = [];
end

spikelist = cmb.WhichClusters(DATA.toplevel);
if spikelist == -1 %no clusters selected
    it = findobj(DATA.toplevel,'Tag','UseCluster1');
    set(it,'value',1);
    spikelist = cmb.WhichClusters(DATA.toplevel);
end

nprobes = length(DATA.probelist);
if length(DATA.probelist) > 1
    DATA.state.includeprobename = 1;
else
    DATA.state.includeprobename = 0;
end
lfpfile = strrep(DATA.datafilename,'.mat','.lfp.mat');
if DATA.bysuffix
    LFP = [];
elseif listonly == 0  && lfponly >= 0 && exist(lfpfile,'file')
    fprintf('Loading %s ... ',lfpfile);
    tic;
    load(lfpfile);
    fprintf(' Took %.1f\n',toc);
end
probedone = DATA.probe;

chk = 0;

if oneprobe | lfponly > 0
    probelist = probedone;
    oneprobe = 1;
else
    probelist = DATA.probelist;
end

d = dir(DATA.datadir);
if nprobes > 1
    str = ['.p' num2str(probedone) 'c1.'];
else
    str = ['.c1.'];
end
nf = 0;
for j = 1:length(d)
    if strfind(d(j).name,str)
        nf = nf+1;
        filenames{nf} = d(j).name;
        xfile(nf) = 1;
    end
end

DATA.state.redoautocut = 1;
oldstate.autoplotnewprobe = DATA.state.autoplotnewprobe;
DATA.state.autoplotnewprobe = 0;


%recombine only expts selected, unless 'All' is selected
expts = get(DATA.clst,'value');
if expts(1) == 1
    expts = 2:length(DATA.exptypelist);
else
    xfile = zeros(size(expts));
end

if nf > 0
    for j = expts
        %    DATA = cmb.ListSubExpts(DATA,j);
        [dp, poutname] = fileparts(cmb.CombinedName(DATA,j,DATA.spikelist(1),'probe',probedone));
        id = strmatch(poutname,filenames);
        if length(id)
            xfile(id) = 0;
        end
    end
    xid = find(xfile); %files with no match
else
    xfile  = zeros(size(DATA.exptypelist));
end



if ~isempty(recombinenames)
    expts = [];
    for j = 1:length(recombinenames)
        s = recombinenames{j};
%first check explabel - rds.dxXce  format        
        id = find(strncmp(s,DATA.explabel(2:end),length(s)));
        if isempty(id)
%now check rds.AC format. N.B. .AC. should do all matching expts            
            if s(end) ~= '.'
                s = [s '.'];
            end
            if s(1) ~= '.' 
                s = ['.' s];
                id = find(CellToMat(strfind(DATA.combinenames(2:end),s)));
                if isempty(id) %in case combine name missing file prefix because state.online set
                    s(1) = '/';
                end
                id = find(CellToMat(strfind(DATA.combinenames(2:end),s)));                
            end
            id = find(CellToMat(strfind(DATA.combinenames(2:end),s)));
        end
        if isempty(id)
            a = split(recombinenames{j},'\.');
            if length(a) < 2
                DATA = AddError(DATA,'-show','RecombineName %s does not have stimulus and type ',recombinenames{j});
            else
            id = find(strcmp(a{2},DATA.expnames));
            if ~isempty(id)
                names = unique(DATA.expstrs(id));
                id = [];
                for k = 1:length(names)
                    name = [a{1} '.' names{k}];
                    b = find(strncmp(name,DATA.explabel(2:end),length(name)));
                    id = [id b(:)'];
                end
            end
            end
        end
        if isempty(id)
            DATA = AddError(DATA,'-show','No Expt %s in %s',recombinenames{j},DATA.name);
        end
        expts = [expts 1+id(:)'];
    end
    nex = length(expts);
elseif onetype == 1
    id = strmatch('orXob',DATA.exptypelist(2:end));
    if isempty(id)
        return;
    else
        onetype = id(1)+1;
    end
    nex = 1;
else
    nex = length(expts)+sum(xfile);
end

if docells %just use P1 to determine expts,
    probelist = 1;
end

for p = 1:length(probelist)
    if listonly == 0 && oneprobe == 0
        if probelist(p) ~= DATA.probe
            DATA = cmb.SetProbe(DATA, probelist(p));
        end
    end
    
    for jx = nex:-1:1
        j = expts(jx);
        if j <= length(DATA.exptypelist)
            set(DATA.clst,'value',j);
            DATA = cmb.ListSubExpts(DATA,j);
            poutname = cmb.CombinedName(DATA,j,DATA.spikelist(1),'probe',probedone);
            donename = cmb.CombinedName(DATA,j,DATA.spikelist(1),'probe',1);
            outname = cmb.CombinedName(DATA,j,DATA.spikelist(1));
            set(DATA.saveitem,'string',outname);
            if ~isempty(DATA.allcombineids{j})
                DATA.combineids = DATA.allcombineids{j};
                fprintf('%s',DATA.exptypelist{j})
                fprintf(',%d',DATA.combineids);
                fprintf('\n');
            else
                DATA.combineids = [];
            end
        else
            poutname = [dp '/' filenames{xid(j-length(DATA.exptypelist))}];
            outname = regexprep(poutname,'\.p[0-9]*c[0-9]\.',sprintf('.p%dc1.',probelist(p)));
            donename = regexprep(poutname,'\.p[0-9]*c[0-9]\.',sprintf('.p1c1.',probelist(p)));
            set(DATA.saveitem,'string',outname);
            DATA.combineids = [];
        end
        DATA.outname = outname;
        if isempty(DATA.combineids)
            if exist(poutname,'file') | exist(donename,'file')
                if exist(poutname,'file')
                    load(poutname);
                else
                    load(donename);
                end
                Expt.Header.Name = cmb.BuildName(Expt.Header.Name);
                DATA.Expt = Expt;
                if j > length(DATA.exptypelist)
                    id =strmatch(Expt.Header.expname,DATA.explist,'exact');
                    if length(id) ==1
                        set(DATA.clst,'value',id);
                        DATA = cmb.ListSubExpts(DATA,id);
                    end
                end
                if isfield(Expt.Header,'Combineids')
                    DATA.combineids = Expt.Header.Combineids;
                    DATA.allcombineids{j} = DATA.combineids;
                    fprintf('%s',splitpath(poutname))
                    fprintf(',%d',DATA.combineids);
                    fprintf('\n');
                elseif isfield(Expt.Header,'Combined')
                    DATA.combineids = DATA.expid(Expt.Header.Combined);
                    DATA.allcombineids{j} = DATA.combineids;
                    fprintf('%s',splitpath(poutname))
                    fprintf(',%d',DATA.combineids);
                    fprintf('\n');
                else
                    DATA.combineids = [];
                    fprintf('%s no combineids',splitpath(poutname));
                    if DATA.logfid
                        fprintf(DATA.logfid,'%s no combineids',splitpath(poutname))
                    end
                end
            else
                DATA.combineids = [];
            end
        end
        if length(DATA.combineids)
            DATA.extype = j;
            id = find(ismember(DATA.expid,DATA.combineids));
            if isempty(id) %ids from a saved file, may not match list in GUI
                DATA.expid = DATA.combineids;
                id = find(ismember(DATA.expid,DATA.combineids));
            end
            
            if ~isempty(id)
                if listonly
                    fprintf('Expts ');
                    fprintf('%d ',id);
                    fprintf('of %d\n',length(DATA.expid));
                else
                    nstr = length(get(DATA.elst,'String'));
                    if max(id) > nstr
                        fprintf('%d is longer than expt list (%d)\n',max(id),nstr);
                    else
                        set(DATA.elst,'value',id);
                    end
                    if lfponly > 0
                        DATA.Expt = cmb.CombinePlot(DATA, 0);
                        if j > length(DATA.exptypelist)
                        end
                        cmb.combine('savelfp',DATA,LFP);
                    else
                        [Expt, DATA] = cmb.CombinePlot(DATA, chk);
                        %             Expt.Header.SpkStats = cmb.GetSpkStats(DATA);
                        drawnow;
                        nspk = sum([Expt.Trials.count]);
                        if j > length(DATA.exptypelist)
                            file = strrep(poutname,['.p' num2str(probedone)],['.p' num2str(probelist(p))]);
                            c = '*'
                        else
                            file = cmb.CombinedName(DATA,j,1);
                            c = '';
                        end
                        if savefile && docells == 0
                            save(file,'Expt');
                            fprintf('Saved %d spikes (Expts %s) to %s (%s)\n',nspk,sprintf(' %d',DATA.combineids),file,DATA.user);
                        end
                        if DATA.logfid
                            fprintf(DATA.logfid, '%s,Saved %d spikes (Expts%s%s) to %s (%s) Recombine\n',datestr(now),nspk,c,sprintf(' %d',DATA.combineids),file,DATA.user);
                        end
                        nspk = Expt.Header.nspk;
                        nc = length(Expt.Header.nspk)-1;
                        cl = 1;
                        spikelist = DATA.spikelist;
                        while cl < nc
                            cl = cl+1;
                            if nspk(cl+1) > 10
                                DATA.spikelist = cl;
                                cmb.SetClusterCheck(DATA);
                                [Expt, DATA] = cmb.CombinePlot(DATA, chk);
                                if j > length(DATA.exptypelist)
                                    file = strrep(poutname,['.p' num2str(probedone) 'c1'],['.p' num2str(probelist(p)) 'c' num2str(cl)]);
                                    c = '*'
                                else
                                    file = cmb.CombinedName(DATA,j,cl);
                                    c = '';
                                end
                                if savefile
                                    save(file,'Expt');
                                end
                            end
                        end
                        DATA.spikelist = spikelist;
                        cmb.SetClusterCheck(DATA);
                        if p == 1 && docells
                            out{jx} = cmb.CombineAllCells(a, DATA, cargs{:});
                        end
                        if p ==1  && lfponly == 0 && (exist('LFP') || docells)
                            DATA.Expt = Expt;
                            cmb.combine('savelfp',DATA,LFP);
                        end
                        %                cmb.CombineAll(a,DATA);
                    end
                end
            end
        else
            fprintf('No combine file list for %s. Using Individual Cell Files.\n',DATA.outname);
            if p == 1 && docells
                out{jx} = cmb.CombineAllCells(a, DATA, cargs{:});
            end
        end
    end
    %Calculate spike shape AFTER combining, so that autocuts are
    %shown
    if onetype == 0 && listonly == 0 && DATA.bysuffix == 0
        DATA = cmb.CalcMeanSpike(DATA,1:length(DATA.Expts));
        cmb.SetExptClusters(DATA);
        if oneprobe ==0 && lfponly < 1
            cmb.SaveSpikeShape(DATA,DATA.meanspkfile);
        end
    end
end
if mkall
    PlotAllProbes(fileparts(DATA.datafilename),'sptrig','save');
end
DATA.state.autoplotnewprobe = oldstate.autoplotnewprobe;


set(DATA.toplevel,'UserData',DATA);

