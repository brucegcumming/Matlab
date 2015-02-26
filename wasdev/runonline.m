function [out, n] = runonline(list, varargin)

% runonline(list)
%
% plots online data files via matlab structures, and allows 
% these matlab structures to be saved.
% Expt = runonline('getexpt')  returns an expt struct of the currenlty
% plotted data

global TOPTAG;
if isempty(TOPTAG)
    TOPTAG = 'runonlineTOPlevel';
end

octc = 0;

j = 1;
while j < nargin
    if ischar(varargin{j})
    if strncmpi(varargin{j},'Tag',3)
        TOPTAG = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'function',7)
        OTTF.xfunc = varargin{j+1};
        j = j+1;
    end
end
    j = j+1;
end


if(isempty(findobj('Tag',TOPTAG)))
  if ~exist('list','var') | isempty(list)
    list = './cedtlist.mall';
else
   % TOPTAG = list;
  end
    OTTF.tag.top = TOPTAG;
     OTTF.tag.fig = 'OTTFDataPlot';
     OTTF.tag.extra = 'OTTFExtraPlot';
     OTTF.tag.figb = 'OTTFDataBPlot';
     OTTF.tag.pop = 'OTTFPopPlot';
     OTTF.tag.select = 'OTTFSelections';
     OTTF.listtag = 'OTTFlist';
     OTTF.savefile = [list '/test.mat'];
     out = InitOTTF(list, OTTF, varargin{:})
else
     top = findobj('Tag',TOPTAG);
     OTTF = get(top,'UserData');
 end
 
 if strcmpi(list,'getstate')
     out = OTTF;
     it = findobj('Tag',OTTF.listtag);
     n = get(it, 'value');
     return;    
elseif strmatch('store',list)
     top = findobj('Tag',TOPTAG);
     set(top,'UserData',varargin{1});
 elseif strmatch('Save',list)
SaveMatFile(OTTF);
elseif strmatch('Refresh',list)
     OTTF.dir = get(findobj('Tag','Prefix'),'String');
     id = findstr(OTTF.dir,'/');

%
%
% HORRIBLE kludge. Matlab cannot execute "system" commands if the current
% path is a network drive (obvious really), so need to temporarily change
% directory while the command is executed....
     
     if ispc
         olddir = pwd;
        cd C:
     end
     if isempty(id)
         if octc
             if ~exist(OTTF.dir)
                 if ispc
                     system(['rsh lsr-mat3 mkdir ' OTTF.dir]);
                 else
                     system(['mkdir ' OTTF.dir]);
                 end
             end

             rdir = strrep(OTTF.dir,'/bgc/data','/local/data');
             if ispc
                 comnd = ['rsh lsr-octd cp ' rdir '/*rc*' OTTF.dir];
                 comb = ['rsh lsr-mat3 rc2tab ' OTTF.dir '/*rc*'];
             else
                 comnd = ['cp ' rdir '/*rc*' OTTF.dir];
                 comb = ['rc2tab ' OTTF.dir '/*rc*'];
             end
             fprintf('%s\n',comnd);
             system(comnd);
             system(comb);
         else
             if ispc
                 tic;
                 system(['rsh lsr-mat3 rc2tab -all ' OTTF.dir]);
                 toc;
             else
                 system(['rc2tab ' OTTF.dir '/*rc*']);
             end
         end
     else
         odir = OTTF.dir(id(1):end);
%         system(['runtab ' odir]); 
         if octc
             rdir = strrep(OTTF.dir,'/bgc/data','/local/data');
             comnd = ['rsh lsr-octc "rc2tab ' rdir '/*rc*"'];
         else
             if ispc
                 comnd = ['rsh lsr-mat3 rc2tab ' odir '/*rc*'];
             else
                 comnd = ['rc2tab ' odir '/*rc*'];
             end
         end
         fprintf('%s\n',comnd);
         system(comnd);    
     end
     if ispc
        cd(olddir);
     end
     OTTF = relist(OTTF);
     runonline('store',OTTF);
 elseif strmatch('loadone',list)
  LoadOTTFData(OTTF.fstrings{OTTF.id});
 elseif strmatch('getexpt',list)
     out = OTTF.ExptData;
     return;
 elseif strmatch('getstate',list)
     out = OTTF;
     return;
 elseif strmatch(list,'load')
  LoadAllData;
elseif strmatch(list,'Mark')
    OTTF.data.marked(OTTF.id) = 1;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
elseif strmatch(list,'UnMark')
    OTTF.data.marked(OTTF.id) = 0;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
    RePlot(OTTF);
elseif strmatch(list,'PrintMarkList')
    idx = find(OTTF.data.marked > 0);
    fprintf('%s\n',OTTF.fstrings{idx});

elseif strmatch(list,'close')
  exit_ottf(OTTF);
elseif strmatch(list,'popselect')
  OTTF.cntrlbox = setselect(OTTF,OTTF.tag.select);
  if ~isempty(OTTF.cntrlbox) && isfield(OTTF,'xfunc')
     feval(OTTF.xfunc,'setselect');
 end
elseif strmatch(list,'TouchPoint')
    j = varargin{2};
    GetFigure(OTTF.tag.extra);
    fprintf('%s\n',OTTF.fstrings{j});
    if OTTF.plot.labels(j) > 0
        delete(OTTF.plot.labels(j));
        OTTF.plot.labels(j) = 0;
    else
        OTTF.plot.labels(j) = text(OTTF.plot.allx(j),OTTF.plot.ally(j),splitpath(OTTF.fstrings{j}));
    end
    OTTF.id = j;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
    doentry(j,OTTF);
elseif strmatch(list,'popplot')
    GetFigure(OTTF.tag.extra);
    hold off;
 for j = 1:length(OTTF.fstrings)
     if exist(OTTF.fstrings{j},'file')
         load(OTTF.fstrings{j});
         quality = CheckSpike(Expt,'noplot');
         plot(quality(1),quality(3),'o','buttondownfcn',['runonline(''TouchPoint'',gcf,' num2str(j) ')']);
         OTTF.plot.allx(j) = quality(1);
         OTTF.plot.ally(j) = quality(3);
         OTTF.plot.labels(j) = 0;
         plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runonline(''TouchPoint'',gcf,' num2str(j) ')']);
         if(OTTF.state.verbose)
             fprintf('%s: %.2f %.2f\n',OTTF.fstrings{j},quality(1),quality(3));
         end
         hold on;
     end
 end
elseif strmatch(list,'next')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n < length(OTTF.listids);
       n = n+1;
       set(it, 'value',n);
       doentry(OTTF.listids(n),OTTF,'setid');
     end
elseif strmatch(list,'prev')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n > 1
       n = n-1;
       set(it, 'value',n);
       doentry(OTTF.listids(n),OTTF,'setid');
     end
 elseif strmatch(list,'replot')
     OTTF.plot.type = get(findobj('Tag','PlotType'),'value') -1;

     OTTF = RePlot(OTTF);
     set(findobj('Tag',TOPTAG),'UserData',OTTF);
     
 elseif strmatch(list,'update')
     OTTF.prefix = get(findobj('Tag','Prefix'),'String');
     OTTF.plot.legendpos = get(findobj('Tag','LegendPos'),'value') -1;
     OTTF.plot.type = get(findobj('Tag','PlotType'),'value') -1;
     OTTF.state.verbose = get(findobj('Tag','Verbose'),'value');
     OTTF.plot.combine = get(findobj('Tag','Combine'),'value');
     OTTF.plot.reverse = get(findobj('Tag','Reverse'),'value');
     OTTF.plot.psych = get(findobj('Tag','PsychAdd'),'value');
     OTTF.plot.autoplot = get(findobj('Tag','Autoplot'),'value');
     if ~isempty(findobj('Tag','ForceZero'))
         OTTF.plot.forcezero = get(findobj('Tag','ForceZero'),'value');
         OTTF.plot.smoothing = get(findobj('Tag','Smooth'),'value');
         OTTF.plot.sdfw = str2num(get(findobj('Tag','sdfw'),'string'));
         s = get(findobj('Tag','rcsdfw'),'string');
         OTTF.plot.rcsdfw = sscanf(s,'%f');
     
         OTTF.plot.nmin = str2num(get(findobj('Tag','Nmin'),'string'));
         OTTF.savefile = get(findobj('Tag','SaveFile'),'string');
         strs = get(findobj('Tag','SDFType'),'string');
         OTTF.plot.sdftype = strs{get(findobj('Tag','SDFType'),'value')};
     end
     runonline('store',OTTF);
     if OTTF.plot.autoplot
         PlotCurrent;
     end

elseif strmatch(list,'setplot')
     SetPlotTypes;
     if OTTF.plot.autoplot
         PlotCurrent;
     end
elseif strmatch(list,'setentry')
    if OTTF.plot.autoplot
        PlotCurrent;
    end
elseif strmatch(list,'Plot')
    PlotCurrent;
elseif strmatch(list,'relist')
     listtype= get(findobj('Tag','ListType'),'value');
     relist(listtype);
end


function OTTF = RePlot(OTTF)

argon = CheckExpt(OTTF.ExptData, {});
if get(findobj('Tag','LFP'),'value')
    argon = {argon{:} {'lfp'}};
end
if get(findobj('Tag','ShowN'),'value')
    argon = {argon{:} {'ShowN'}};
end

if get(findobj('Tag','Pcolor'),'value')
    argon = {argon{:} {'pcolor'}};
end

if strncmpi('rates',OTTF.plot.onetype,5)
    argon = {argon{:} {'rates'}};
end
if strncmp('cycSDF',OTTF.plot.onetype,6)
    argon = {argon{:} {'sdf'} {'periodic'} {'sdfw'} OTTF.plot.sdfw};
end
if strncmp('SDF',OTTF.plot.onetype,3)
    argon = {argon{:} {'sdf'} {'showsum'} {'sdfw'} OTTF.plot.sdfw};
end
if OTTF.plot.psych
    argon = {argon{:} {'psych'}};
end
if strmatch(OTTF.plot.sdftype,'Box')
    argon = {argon{:}, OTTF.plot.sdftype};
end
GetFigure(OTTF.tag.fig);
Expt = OTTF.ExptData;
if ~isempty(OTTF.Txtdata)
    plottype = OTTF.plot.type;
    if OTTF.plot.type == 1
        NsineRC(OTTF.ExptData,OTTF.Txtdata,'monoc');
    elseif ismember(plottype,[5 6 7]) %%runnsine plots that need sdf
        if ~isfield(OTTF,'nsdfs')
            fprintf('Building SDFs\n');
%            OTTF.nsdfs = runnsine(OTTF.Txtdata,Expt,'mksdf');
OTTF.nsdfs = [];
        end
        if plottype == 5
            res = runnsine(OTTF.Txtdata,Expt,OTTF.nsdfs,'sinamps');
        elseif plottype == 6
            res = runnsine(OTTF.Txtdata,Expt,OTTF.nsdfs,'sfdp');
        end
        if isfield(res,'nsdfs') & isempty(OTTF.nsdfs)
            OTTF.nsdfs = res.nsdfs;
        end
    else
        if isfield(OTTF.Txtdata,'delays')
            NsineRC(Expt,OTTF.Txtdata);
        else
            nsres = NsineRC(Expt,OTTF.Txtdata,[],{});
            OTTF.Txtdata.delays = nsres.delays;
            OTTF.Txtdata.rcid = nsres.rcid;
        end
    end
elseif(Expt.isrc)
%             if strmatch(Expt.Stimvals.et,'dO')
%                 PlotRevCor(Expt,'sdfw',100,argon{:});
 %                PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,argon{:});
 %           else
 if OTTF.plot.type == 2
     PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,'nmin',OTTF.plot.nmin,argon{:},'yval',0,'psych');
 elseif OTTF.plot.type == 3
     PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,'nmin',OTTF.plot.nmin,argon{:},'yval',0);
 elseif OTTF.plot.type == 8
     for j = 1:length(Expt.Trials)
         durs(j) = Expt.Trials(j).End(1)-Expt.Trials(j).Start(1);
     end
     duration = min(durs);
     if duration < 100
         duration = 22000;
     end
     times = [0:10:duration];
     sdf = trigsdf(Expt.Trials,OTTF.plot.rcsdfw,times);
     subplot(1,2,1);
     plot(times,sdf);
     title(sprintf('TF %.2f',Expt.Stimvals.tf));
     for j = 1:length(Expt.Trials)
         Expt.Trials(j).Trigger = Expt.Trials(j).Start(1);
     end
     period = 10000/Expt.Stimvals.tf;
     times = [-period:10:period];
     sdf = trigsdf(Expt.Trials,OTTF.plot.rcsdfw,times,'event');
     subplot(1,2,2);
     hold off;
     plot(times./10,sdf);
     hold on;
     for j = 1:length(Expt.Trials)
         Expt.Trials(j).Trigger = Expt.Trials(j).Start(1);
     end
     period = 10000/Expt.Stimvals.tf;
     times = [0:10:duration];
     sdf = trigsdf(Expt.Trials,OTTF.plot.rcsdfw,times,'period',period);
     plot([1:length(sdf)]./10,sdf,'r');
 else
     PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,'nmin',OTTF.plot.nmin,argon{:});
 end
 %end
elseif(OTTF.plot.reverse)
    OTTF.data.plotdata{OTTF.id} = PlotRates(Expt, Expt.Stimvals.e2, 'type2',Expt.Stimvals.et, 'legendpos',OTTF.plot.legendpos,argon{:});
else
    OTTF.data.plotdata{OTTF.id} = PlotRates(Expt, Expt.Stimvals.et, 'type2',Expt.Stimvals.e2, 'legendpos',OTTF.plot.legendpos,argon{:});
end

function doentry(id, OTTF,varargin)
   
j = 1;
while(j < nargin - 1)
    if strncmpi(varargin{j},'setid',5)
        OTTF.id = 1;
    end
    j = j+1;
end

argon = {};
if get(findobj('Tag','LFP'),'value')
    argon = {argon{:} {'lfp'}};
end
if get(findobj('Tag','ShowN'),'value')
    argon = {argon{:} {'ShowN'}};
end

if get(findobj('Tag','Pcolor'),'value')
    argon = {argon{:} {'pcolor'}};
end

if strncmpi('rates',OTTF.plot.onetype,5)
    argon = {argon{:} {'rates'}};
end
if strncmp('cycSDF',OTTF.plot.onetype,6)
    argon = {argon{:} {'sdf'} {'periodic'} {'sdfw'} OTTF.plot.sdfw};
end
if strncmp('SDF',OTTF.plot.onetype,3)
    argon = {argon{:} {'sdf'} {'showsum'} {'sdfw'} OTTF.plot.sdfw};
end

GetFigure(OTTF.tag.fig);
Txt = [];
np = 1;
if OTTF.plot.combine
    for j = 1:length(id)
        mafile = strrep(OTTF.fstrings{id(j)},'.st','.ma');
        if OTTF.state.verbose
            fprintf('Reading %s\n',mafile);
        end
        exps{np} = ReadOnline(mafile);
        if isfield(exps{np},'txtable');
            tables{np} = dlmread(exps{np}.txtable);
            Txt = [Txt; tables{np}];
        end
        np = np+1;
    end
    if OTTF.state.verbose
        fprintf('Combining files\n');
    end
    Expt = AddExpts(exps{:});
    
    argon = CheckExpt(Expt,argon);
    hold off;
    OTTF.ExptData = Expt;
    if OTTF.state.verbose
        fprintf('Plotting\n');
    end
    if ~isempty(Txt)
        im = tab2im(Txt,Expt);
        OTTF.Txtdata = im;
        OTTF.nsdfs = [];
    else
        OTTF.Txtdata = [];
    end
    OTTF = RePlot(OTTF);
else
    alltitle = [];
    for j = 1:length(id)
        mafile = strrep(OTTF.fstrings{id(j)},'.st','.ma');
        if OTTF.state.verbose
            fprintf('Reading %s\n',mafile);
        end
        Expt = ReadOnline(mafile);
        argon = CheckExpt(Expt,argon);
        if np > 1
            hold on;
        else
            hold off;
        end
        if OTTF.state.verbose
            fprintf('Plotting\n');
        end
        if(Expt.isrc)
            PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,'nmin',OTTF.plot.nmin,argon{:});
        else
            OTTF.data.plotdata{id(j)} = PlotRates(Expt, Expt.Stimvals.et, 'type2',Expt.Stimvals.e2, 'legendpos',OTTF.plot.legendpos,argon{:});
            alltitle = [alltitle sprintf('%s at %.0f\n',OTTF.data.plotdata{id(j)}.title,Expt.Trials(1).Start)];
        end
        np = np+1;
    end
    OTTF.ExptData = Expt;
    if np < 4
        title(alltitle)
    else
        title(sprintf('%s',OTTF.data.plotdata{id(j)}.title));
    end
end
runonline('store',OTTF);

 function SaveMatFile(OTTF)
 
  Expt = OTTF.ExptData;
  save(OTTF.savefile,'Expt');

 function args = CheckExpt(Expt,args)

 if Expt.Stimvals.st == 2 || Expt.Stimvals.st == 15  %% RDS or RLS
     args = {args{:} 'Uncorr'};
 end
 if strcmp(Expt.Stimvals.et,'sf') | strcmp(Expt.Stimvals.et,'tf')
     args = {args{:} 'LogX'};
 end
 if strmatch(Expt.Stimvals.et,{'dp', 'dO', 'dx'}) & ~ Expt.isrc
     args = {args{:} 'AddMonoc'};
 end

function PlotCurrent()
    [OTTF, n] = runonline('getstate');

    doentry(OTTF.listids(n),OTTF,'setid');

function SetPlotTypes()
%SetPlotTypes looks at the values of the menus for selecting
%plot types, and records them in OTTF for use elsewhere

OTTF = runonline('getstate');

it = findobj('Tag','plottype');
type = get(it, 'value');
str = get(it, 'String');
OTTF.plot.onetype = strrep(str(type,:),' ',''); %The string shown on the menu.

it = findobj('Tag','extraplot');
if ~isempty(it)
    type = get(it, 'value');
    str = get(it, 'String');
    OTTF.plot.extraplot = strrep(str(type,:),' ',''); %The string shown on the menu.
end
doentry(OTTF.id,OTTF);
runonline('store',OTTF);


function LoadAllData()
    global OTTF;

    for j = 1:length(OTTF.fstrings)
        fprintf('%s',OTTF.fstrings{j});
        LoadOTTFData(OTTF.fstrings{j});
    end

function CloseTag(tag)
    it = findobj('Tag',tag);
    if ~isempty(it)
        close(it);
    end


function exit_ottf(OTTF)
    names = fieldnames(OTTF.tag);
    for j = 1:length(names)
        eval(['CloseTag(OTTF.tag.' names{j} ');']);
    end
    return;


function OTTF = relist(OTTF)


    lso = dir(OTTF.dir);
    sts = regexp({lso.name},'.st[0-9]*$');
    OTTF.strings = {};
    for j = 1:length({lso.name})
        if ~isempty(sts{j})
            if isfield(OTTF,'strings') & ~isempty(OTTF.strings)
                OTTF.strings = {OTTF.strings{:} lso(j).name};
            else
                OTTF.strings{1} = lso(j).name;
            end
        end
    end

    for j = 1:length(OTTF.strings)
        OTTF.fstrings{j} = [OTTF.dir '/' OTTF.strings{j}];
        e = ReadOnline(OTTF.fstrings{j},'Header');
        OTTF.strings{j} = sprintf('%s %.1f %.1f',OTTF.strings{j},e.Stimvals.ns,e.Stimvals.nr);
        if isfield(e.Stimvals,'pt') & e.Stimvals.pt > 1
            OTTF.strings{j} = [OTTF.strings{j} sprintf(' %dRC',e.Stimvals.nf+1)];
        end        
        if isfield(e.Stimvals,'sl') & e.Stimvals.sl == 1
            OTTF.strings{j} = [OTTF.strings{j} 'sl1'];
        end
        OTTF.listids(j) = j;
    end

    lst = findobj('Tag',OTTF.listtag)
    set(lst, 'String',OTTF.strings,'value',1);

function OTTF = InitOTTF(odir, OTTF, varargin)
global bgcfileprefix;

if ~exist(odir,'dir')
    fprintf('No Dir %s',odir);
    return;
end


if ispc
    tic;
    olddir = pwd;
    cd C: 
    comd = ['rsh lsr-mat3 rc2tab -all ' odir]
    system(comd);
    cd(olddir);
    toc;
else
    system(['rc2tab ' odir '/*rc*']);
end

lso = dir(odir);
sts = regexp({lso.name},'.st[0-9]*$');
rcs =  regexp({lso.name},'.rc[0-9]*$');

if isfield(OTTF,'strings')
    OTTF = rmfield(OTTF,'strings');
end

%            system(['runtab ' lso(j).name]);




OTTF.dir = odir;
OTTF = relist(OTTF);


if ~isempty(bgcfileprefix)
    OTTF.prefix = bgcfileprefix;
else
    OTTF.prefix = '';
end

OTTF.listids = 1:length(OTTF.fstrings);
OTTF.state.showonset= 0;
OTTF.state.minsacs = 1;
OTTF.state.minspikes = 1;
OTTF.state.mintrials = 1;
OTTF.state.minsacsize = 0.1;
OTTF.state.maxsacsize = 1.0;
OTTF.state.addedspikes = 0;
OTTF.state.invertspikes = 0;
OTTF.state.rebuild = 0;
OTTF.state.verbose = 0;
OTTF.state.includeratio = 0.5;
OTTF.state.dirw = 0.5;
OTTF.plot.onetype = 'Saccade';
OTTF.plot.type = 0;
OTTF.plot.nmin = 100;
OTTF.plot.tftype = 1;
OTTF.plot.sizetype = 1;
OTTF.plot.poptype = 1;
OTTF.plot.sdfw = 150;
OTTF.plot.rcsdfw = 100;
OTTF.plot.sdftype = 'Exp';
OTTF.plot.sdftypes = {'Exp' 'Box' 'Gauss' 'HalfG'};
OTTF.plot.legendtype = 1;
OTTF.plot.sdfstart = -5000;
OTTF.plot.sdfend = 5000;
OTTF.plot.sdfstep = 20;
OTTF.plot.legendpos = 2;
OTTF.plot.extraplot = 'None';
OTTF.plot.forcezero = 0;
OTTF.plot.autoplot = 0;
OTTF.plot.smoothing = 0;
OTTF.plot.combine = 1;
OTTF.plot.reverse = 0;
OTTF.plot.psych = 0;
j = 1;
while(j < nargin-1)
    if strncmpi(varargin{j},'prefix',5)
        j = j+1
        OTTF.prefix = varargin{j};
    elseif strncmpi(varargin{j},'verbose',5)
        j = j+1
        OTTF.state.verbose = 1;
    end
    j = j+1;
end


if ~isempty(OTTF.prefix)
    for j = 1:length(OTTF.fstrings)
        OTTF.fstrings{j} = [OTTF.prefix OTTF.fstrings{j}];
    end
end

for j = 1:length(OTTF.fstrings)
    OTTF.data.marked(j) = 0;
end

idx = findstr(odir,'/');
if isempty(idx)
  idx = findstr(odir,'\\');
end
if isempty(idx)
  listname = odir;
else
  listname = odir(idx(end)+1:end);
end

  scrsz = get(0,'Screensize');
  boxw = 300;
  cntrl_box = figure('Position', [100 scrsz(4)-340 boxw 310],...
       'NumberTitle', 'off', 'Tag',OTTF.tag.top,'Name',sprintf('Online Data %s (%d)',listname,length(OTTF.fstrings)));
  lst = uicontrol(gcf, 'Style','listbox','String',OTTF.strings,...
      'Max',10,...
		'Callback', ' runonline(''setentry'')','Tag',OTTF.listtag,...
		'Position',[10 10 boxw-20 150]);

  SPACE = 3;
  VSPACE = 10;
  cw = 10;

  bp(1) = SPACE; bp(2) = 160; bp(3) = 40; bp(4) = 22;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runonline(''next'')',...
'String', '>>', 'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runonline(''prev'')',...
'String', '<<', 'Position', bp);



bp(1) = bp(1) + bp(3)+SPACE;
bp(3) = 40;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runonline(''Refresh'')',...
'String', 'Refresh', 'Position', bp);
bp(1) = bp(1) + bp(3)+SPACE;
bp(3) = 40;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runonline(''Plot'')',...
'String', 'Plot', 'Position', bp);
bp(1) = bp(1) + bp(3)+SPACE;
bp(3) = 40;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runonline(''Save'')',...
'String', 'Save', 'Position', bp);
bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 70;
uicontrol(gcf,'style','pop','string','Counts|Rates|SDF|cycSDF', ...
    'Callback', ' runonline(''setplot'')', 'Tag','plottype',...
    'position',bp,'value',1);


bp(1) = SPACE;
bp(2) = bp(2)+bp(4)+VSPACE;
bp(3) = boxw-SPACE;
uicontrol(gcf,'Style', 'edit', 'Callback', 'runonline(''Refresh'')',...
'String', OTTF.dir, 'Tag', 'Prefix','Position', bp);


    bp(1) = SPACE; bp(3) = 25; bp(2) = bp(2) + bp(4)+VSPACE; bp(4) = 22;
bp(3) = 50;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runonline(''update'')',...
'String', 'Auto', 'Tag', 'Autoplot', 'Position', bp,'value',OTTF.plot.autoplot);


  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 6 * cw;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runonline(''update'')',...
'String', 'ShowN', 'Tag', 'ShowN', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 5*cw;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Pcolor', 'Tag', 'Pcolor', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 7 * cw;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runonline(''update'')',...
'String', 'Verbose', 'Tag', 'Verbose', 'Position', bp,'value',OTTF.state.verbose);

bp(1) = SPACE;
bp(2) = bp(2)+bp(4)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runonline(''update'')',...
'String', 'Combine', 'Tag', 'Combine', 'Position', bp,'value',OTTF.plot.combine);
bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runonline(''update'')',...
'String', 'Psych', 'Tag', 'PsychAdd', 'Position', bp,'value',OTTF.plot.psych);

  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw * 6;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runonline(''update'')',...
'String', 'Reverse', 'Tag', 'Reverse', 'Position', bp,'value',OTTF.plot.reverse);

bp(1) = SPACE;
bp(2) = bp(2)+bp(4)+SPACE;
  uicontrol(gcf,'Style', 'text','String','Legend','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw * 5;
  uicontrol(gcf,'style','pop','string','Auto|UR|UL|LL|LR|Off|None', ...
		    'Callback', ' runonline(''update'')', 'Tag','LegendPos',...
		    'position',bp);
  
  bp(1) = bp(1) +bp(3) + SPACE;
  uicontrol(gcf,'Style', 'text','String','Type','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw *7;
  uicontrol(gcf,'style','pop','string','Default|Monocs|PsyRC|RC0|NSsina|NSdpamp|NSdpsf|Psych|Periodic|New', ...
		    'Callback', ' runonline(''replot'')', 'Tag','PlotType',...
		    'position',bp);

        x = 10; y = 180;

  hm = uimenu(gcf,'Label','Mark');
  uimenu(hm,'Label','&Mark','Callback','runonline(''Mark'')');
  uimenu(hm,'Label','UnMark','Callback','runonline(''UnMark'')');
  uimenu(hm,'Label','List','Callback','runonline(''PrintMark'')');
  uimenu(hm,'Label','UnMark All','Callback','runonline(''UnMarkAll'')');
  uimenu(hm,'Label','Load One','Callback',' runonline(''loadone'')');
  uimenu(hm,'Label','Load All','Callback',' runonline(''load'')');
  uimenu(hm,'Label','Close','Callback',' runonline(''close'')');
  hm = uimenu(gcf,'Label','Options');
  uimenu(hm,'Label','Options','Callback',' runonline(''popselect'')');
  set(gcf,'Menubar','none');
  hm = uimenu(gcf,'Label','&Next','Callback','runonline(''next'');');

OTTF.figid.mainplot = figure('Tag',OTTF.tag.fig,'Position',[300 scrsz(4)-470 512 ...
		    400]);
OTTF.id = 1;
set(cntrl_box,'UserData',OTTF);

function data = Condense(res)


data.sf = res.x;
data.n = res.n;
data.eye = res.y;
data.count = res.means;
data.sd = res.sd;
data.name = res.name;
data.extras = res.extras;
data.width = GetEVal(res.Data,'wi');
data.height = GetEVal(res.Data,'hi');
data.ori = GetEVal(res.Data,'or');
data.x = GetEVal(res.Data,'xo');
data.y = GetEVal(res.Data,'yo');


fz = GetEval(res.Data,'fz');
if isnan(fz)
    fz = 72.24;
end
dur = mean([res.Data.Trials.End] - [res.Data.Trials.Start]);
data.duration = round(dur * fz/10000) * 1000/fz;


function OTTF = GetSimple(OTTF, varargin)


[a, b, names] = textread('../tf/simplelist','%f %f %s');

for j = 1:length(OTTF.fstrings)
    name = splitpath(OTTF.fstrings{j});
    OTTF.data.sxcx(j) = NaN;
    for k = 1:length(names)
        if strncmp(names{k},name,11)
            OTTF.data.sxcx(j) = b(k);
        end
    end
end
runonline('store',OTTF);


function cntrl_box = setselect(OTTF, tag)

cntrl_box = [];

if isempty(findobj('Tag',tag))
    
    if isempty(OTTF.plot.forcezero)
        OTTF.plot.forcezero = 0;
    end
    if isempty(OTTF.plot.autoplot)
        OTTF.plot.autoplot = 0;
    end
    if isempty(OTTF.plot.smoothing)
        OTTF.plot.smoothing = 0;
    end

  SPACE = 5;
  VSPACE = 2;
  h = 220;
  scrsz = get(0,'Screensize');
 wsc = scrsz(3) /1000;
  cntrl_box = figure('Position', [200 scrsz(4)-(h+30)*wsc 280*wsc h*wsc], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag',tag,'Name','Section Criteria');

   cw = 10*wsc;
   ch = 18;
   bh = 18*wsc;
   bp(1) = SPACE;
   bp(2) = ch+VSPACE;
   bp(3) = cw*7;
   bp(4) = ch;
  
  uicontrol(gcf,'Style', 'CheckBox','String','Force Zero','Position', bp,...
      'Tag','ForceZero','Callback','runonline(''update'')','value',OTTF.plot.forcezero);
  bp(1) = bp(1) + bp(3) + SPACE;
  bp(3) = cw*5;
    uicontrol(gcf,'Style', 'CheckBox','String','Smooth','Position', bp,...
      'Tag','Smooth','Callback','runonline(''update'')','value',OTTF.plot.smoothing);

  bp(1) = SPACE;
  bp(2) = bp(2) + ch + VSPACE;
  bp(3) = cw*4;
  uicontrol(gcf,'Style', 'text','String','sdfw','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' runonline(''update'')',...
    'String', sprintf('%d',OTTF.plot.sdfw), 'Tag', 'sdfw','Position', bp);

  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*6;
  uicontrol(gcf,'Style', 'text','String','RCsdfw','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 3 * cw;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' runonline(''update'')',...
    'String', sprintf('%d',OTTF.plot.rcsdfw), 'Tag', 'rcsdfw','Position', bp);

  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*6;
  uicontrol(gcf,'Style', 'text','String','Nmin','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 3 * cw;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' runonline(''update'')',...
    'String', sprintf('%d',OTTF.plot.nmin), 'Tag', 'Nmin','Position', bp);

bp(1) = SPACE;
  bp(2) = bp(2) + ch + VSPACE;
  bp(3) = cw*3;
  uicontrol(gcf,'Style', 'text','String','SDF:','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 4 * cw;
  
  val = strmatch(OTTF.plot.sdftype,OTTF.plot.sdftypes);
  uicontrol(gcf,'Style', 'pop', 'Callback', ' runonline(''update'')',...
    'String', OTTF.plot.sdftypes, 'Tag', 'SDFType','Position', bp,'value',val);

  
  bp(1) = SPACE;
  bp(2) = bp(2) + ch + VSPACE;
  bp(3) = cw*6;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runonline(''Save'')',...
'String', 'Save As', 'Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 20 * cw;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' runonline(''update'')',...
    'String', OTTF.savefile, 'Tag', 'SaveFile','Position', bp);

OTTF.wsc = wsc;
  OTTF.cntrlbox = cntrl_box;
  OTTF.selectlast = bp(2);
  runonline('store',OTTF);
end
