function [out, n] = runlist(list, varargin)

% Basic list -running GUI.
% Two other functions required separately are
%
% doentry (id, OTTF) does requested action (plot,fit, etc) on selected file
% setselect(OTTF) options setting panel.T
%
% all the state variables are stored in the top figure, in OTTF

global TOPTAG;
if isempty(TOPTAG)
    TOPTAG = 'runlistTOPlevel';
    tprefix = 'OTTF';
else
    tprefix = TOPTAG;
end


j = 1;
while j < nargin
    if ischar(varargin{j})
    if strncmpi(varargin{j},'Tagpref',6)
        tprefix = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'Tag',3)
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
  OTTF.tag.fig = strrep('OTTFDataPlot','OTTF',tprefix);
  OTTF.tag.extra = 'OTTFExtraPlot';
  OTTF.tag.figb = 'OTTFDataBPlot';
  OTTF.tag.pop = 'OTTFPopPlot';
  OTTF.tag.select = 'OTTFSelections';
  OTTF.tag.data = 'OTTFdata';
  OTTF.listtag = 'OTTFlist';
  out = InitOTTF(list, OTTF, varargin{:});
  OTTF = out;
else
    top = findobj('Tag',TOPTAG);
    if strmatch('store',list)
        set(top,'UserData',varargin{1});
        return;
    end
    OTTF = get(top,'UserData');
end

j = 1;
while j <= length(varargin)
    if ischar(varargin{j})
        if strncmpi(varargin{j},'load',4)
            OTTF = LoadState(OTTF);
            set(OTTF.toplevel,'UserData',OTTF);
            out = OTTF;
        end
    end
    j = j+1;
end


 if strcmpi(list,'getstate')
     out = OTTF;
     it = findobj(OTTF.toplevel,'Tag',OTTF.listtag);
     n = get(it, 'value');
     return;    
elseif strmatch('store',list)
     top = findobj('Tag',TOPTAG);
     set(top,'UserData',varargin{1});
elseif strmatch('loadone',list)
  LoadOTTFData(OTTF,OTTF.id);
elseif strmatch(list,'load')
  LoadAllData(OTTF);
elseif strmatch(list,'Mark')
    OTTF.data.marked(OTTF.id) = 1;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
 elseif strmatch(list,'NewPrefix')
     OTTF.prefix = get(findobj('Tag','Prefix'),'String');
     for j =1:length(OTTF.strings)
         OTTF.fstrings{j} =  [OTTF.prefix OTTF.strings{j}];
     end
    set(findobj('Tag',TOPTAG),'UserData',OTTF);     
elseif strmatch(list,'UnMarkAll')
    OTTF.data.marked(1:end) = 0;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
    RePlot(OTTF);
elseif strmatch(list,'UnMark')
    OTTF.data.marked(OTTF.id) = 0;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
    RePlot(OTTF);
elseif strmatch(list,'PrintMarkList')
    idx = find(OTTF.data.marked > 0);
    if isempty(idx) && isfield(OTTF.plot,'selected')
        idx = OTTF.plot.selected;
    end
    fprintf('%s\n',OTTF.fstrings{idx});
    if isfield(OTTF,'xfunc')
        OTTF = feval(OTTF.xfunc,'DATA',OTTF,varargin{:},list);
    end
elseif strmatch(list,{'All' 'More'})
    if strmatch(list,'All')
        startid = 1;    
    else
        it = findobj(gcf, 'Tag',OTTF.listtag);
        n = get(it, 'value');
        startid = n;
    end
    if isfield(OTTF,'xfunc')
        OTTF = feval(OTTF.xfunc,'DATA',OTTF,varargin{:},list);
        OTTF.test = 1;
    elseif strmatch(OTTF.plot.xop,'MakeSF')
        OTTF = GetSimple(OTTF);
        GetFigure(OTTF.tag.fig);
        clear SFData;
        for j = startid:length(OTTF.fstrings)
            res = PlotExpt(OTTF.fstrings{j});
            drawnow;
            SFData(j) = Condense(res);
        end
        for j = 1:length(OTTF.fstrings)
            SFData(j).complex = OTTF.data.sxcx(j);
        end
        save('SFData.mat','SFData');
    else
        for j = startid:length(OTTF.fstrings)
            OTTF = doentry(j,OTTF);
            drawnow;
        end
    end
elseif strmatch(list,'close')
  exit_ottf(OTTF);
elseif strmatch(list,'stripprefix')
for j = 1:length(OTTF.fstrings)
    OTTF.fstrings{j} = strrep(OTTF.fstrings{j},'C:\','/');    
end
     runlist('store',OTTF);

     
 elseif strmatch(list,'popselect')
  OTTF.cntrlbox = setselect(OTTF,OTTF.tag.select);
  if ~isempty(OTTF.cntrlbox) && isfield(OTTF,'xfunc')
     feval(OTTF.xfunc,'setselect');
 end
elseif strmatch(list,'TouchPoint')
    j = varargin{2};
    GetFigure(OTTF.tag.pop);
    fprintf('%s\n',OTTF.fstrings{j});
    if OTTF.plot.labels(j) > 0
        delete(OTTF.plot.labels(j));
        OTTF.plot.labels(j) = 0;
    else
        if ~strmatch(OTTF.plot.onetype,'Replay')
        OTTF.plot.labels(j) = text(OTTF.plot.allx(j),OTTF.plot.ally(j),splitpath(OTTF.fstrings{j}));
        end
    end
    OTTF.id = j;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
    doentry(j,OTTF);
elseif strmatch(list,'popplot')
    GetFigure(OTTF.tag.pop);
    hold off;
    strs = get(findobj('Tag','PopPlot'),'string');
    val =  get(findobj('Tag','PopPlot'),'value');
    type = strrep(strs(val,:),' ','');
    OTTF.plot.poptype = type;
    if strcmp(type,'Isolation')
        for j = 1:length(OTTF.fstrings)
            if exist(OTTF.fstrings{j},'file')
                load(OTTF.fstrings{j});
                quality = CheckSpike(Expt,'noplot');
                plot(quality(1),quality(3),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                OTTF.plot.allx(j) = quality(1);
                OTTF.plot.ally(j) = quality(3);
                OTTF.plot.labels(j) = 0;
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                if(OTTF.state.verbose)
                    fprintf('%s: %.2f %.2f\n',OTTF.fstrings{j},quality(1),quality(3));
                end
                hold on;
            end
        end
    elseif strmatch(type,{'SpikeW'})
        for j = 1:length(OTTF.fstrings)
            if exist(OTTF.fstrings{j},'file')
                load(OTTF.fstrings{j});
                ds = SpikeShape(Expt.Spike);
                OTTF.plot.allx(j) = ds.width;
                OTTF.plot.ally(j) = ds.dprime;
                OTTF.plot.labels(j) = 0;
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                if(OTTF.state.verbose)
                    fprintf('%s: %.2f %.2f\n',OTTF.fstrings{j},quality(1),quality(3));
                end
                hold on;
            end
        end
    elseif strmatch(type,{'Psych'})
        if ~isfield(OTTF.plot,'minpsychtrials')
            OTTF.plot.minpsychtrials = 100;
        end
        ratios = [];
        for j = 1:length(OTTF.popdata)
            if ~isempty(OTTF.popdata{j}) & OTTF.popdata{j}.ntrials > OTTF.plot.minpsychtrials
                OTTF.plot.allx(j) = OTTF.popdata{j}.psf(1);
                OTTF.plot.ally(j) = OTTF.popdata{j}.psf(2);
                ratios(j) = OTTF.popdata{j}.psf(2)/OTTF.popdata{j}.psf(1);
                OTTF.plot.labels(j) = 0;
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                hold on;
            else
                ratios(j) = NaN;
            end
        end
        MarkIdentity(gca);
        title(sprintf('Geometric mean %.3f\n',exp(mean(log(ratios(find(~isnan(ratios))))))));
    elseif strmatch(type,{'Psychb'})
        ratios = [];
        for j = 1:length(OTTF.popdata)
            if ~isempty(OTTF.popdata{j})
                ratios(j) = OTTF.popdata{j}.psf(2)/OTTF.popdata{j}.psf(1);
                OTTF.plot.allx(j) = OTTF.popdata{j}.ntrials;
                OTTF.plot.ally(j) = ratios(j);
                OTTF.plot.labels(j) = 0;
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                hold on;
            else
                ratios(j) = NaN;
            end
        end
        set(gca,'Yscale','log');
        title(sprintf('Geometric mean %.3f\n',exp(mean(log(ratios(find(~isnan(ratios))))))));
    elseif strmatch(type,{'New'})
        ratios = [];
        for j = 1:length(OTTF.popdata)
            if ~isempty(OTTF.popdata{j})
                ratios(j) = OTTF.popdata{j}.psf(2)/OTTF.popdata{j}.psf(1);
                OTTF.plot.allx(j) = OTTF.popdata{j}.mean(1);
                OTTF.plot.ally(j) = mean(2);
                OTTF.plot.labels(j) = 0;
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                hold on;
            else
                ratios(j) = NaN;
            end
        end
        set(gca,'Yscale','log');
        title(sprintf('Geometric mean %.3f\n',exp(mean(log(ratios(find(~isnan(ratios))))))));
    elseif strmatch(type,{'ConsecLFP'})
        ratios = [];
        for j = 1:length(OTTF.popdata)
            if ~isempty(OTTF.popdata{j})
                ratios(j) = OTTF.popdata{j}.psf(2)/OTTF.popdata{j}.psf(1);
                OTTF.plot.allx(j) = OTTF.popdata{j}.ntrials;
                OTTF.plot.ally(j) = ratios(j);
                OTTF.plot.labels(j) = 0;
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                hold on;
            else
                ratios(j) = NaN;
            end
        end
        set(gca,'Yscale','log');
        title(sprintf('Geometric mean %.3f\n',exp(mean(log(ratios(find(~isnan(ratios))))))));
    elseif strmatch(type,{'ConsecMeans'})
        ratios = [];
        for j = 1:length(OTTF.popdata)
            if ~isempty(OTTF.popdata{j})
                ratios(j) = OTTF.popdata{j}.psf(2)/OTTF.popdata{j}.psf(1);
                OTTF.plot.allx(j) = OTTF.popdata{j}.ntrials;
                OTTF.plot.ally(j) = ratios(j);
                OTTF.plot.labels(j) = 0;
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                hold on;
            else
                ratios(j) = NaN;
            end
        end
        set(gca,'Yscale','log');
        title(sprintf('Geometric mean %.3f\n',exp(mean(log(ratios(find(~isnan(ratios))))))));
    elseif strmatch(type,{'LFPCoh', 'LFPnspk', 'LFPCohF', 'LFPstim','LFPCS'})
        subtype = strmatch(type,{'LFPCoh', 'LFPnspk', 'LFPCohF', 'LFPstim','LFPCS'});
        subtype = subtype(1);
        for j = 1:2
            avg(j).alltfs = [];
            avg(j).alln = [];
            avg(j).resp = [];
            avg(j).rate = [];
        end
        for j = 1:length(OTTF.fstrings)
            if isfield(OTTF.popdata{j},'depth') & isfield(OTTF.popdata{j},'lfp') & OTTF.popdata{j}.lfp.nspk > 100
                if ismember(subtype, [1 4])
                    OTTF.plot.allx(j) = OTTF.popdata{j}.depth;
                elseif subtype == 2
                    OTTF.plot.allx(j) = OTTF.popdata{j}.lfp.nspk;
                elseif subtype == 3
                    OTTF.plot.allx(j) = OTTF.popdata{j}.lfp.coh_peakf;
                elseif subtype == 5
                    OTTF.plot.allx(j) = OTTF.popdata{j}.lfp.conc;
                else
                    OTTF.plot.allx(j) = OTTF.popdata{j}.lfp.stimcoh;
                end
                if ismember(subtype,[1 2 3])
                    OTTF.plot.ally(j) = OTTF.popdata{j}.lfp.conc;
                    if OTTF.plot.ally(j) > 0.5
                        av = 1;
                    else
                        av = 2;
                    end
                else
                    OTTF.plot.ally(j) = OTTF.popdata{j}.lfp.stimcoh;
                    if OTTF.plot.ally(j) > 1.4
                        av = 1;
                    else
                        av = 2;
                    end
                end

                oldtfs = avg(av).alltfs;
                oldn = avg(av).alln;
                oldresp = avg(av).resp;
                oldrate = avg(av).rate;
                avg(av).alltfs = sort(unique([avg(av).alltfs OTTF.popdata{j}.lfp.tfs]));
                id = find(ismember(avg(av).alltfs,OTTF.popdata{j}.lfp.tfs));
                oid = find(ismember(avg(av).alltfs,oldtfs));
                if isempty(oldn)
                    avg(av).alln = ones(size(avg(av).alltfs));
                    avg(av).resp = OTTF.popdata{j}.lfp.tfgamma;
                    avg(av).rate = OTTF.popdata{j}.lfp.rates./max(OTTF.popdata{j}.lfp.rates);
                else
                    avg(av).alln(oid) = oldn;
                    avg(av).resp = zeros(size(avg(av).alltfs));
                    avg(av).rate = zeros(size(avg(av).alltfs));
                    avg(av).alln = zeros(size(avg(av).alltfs));
                    avg(av).alln(oid) = oldn;
                    avg(av).resp(oid) = oldresp;
                    avg(av).resp(id) = avg(av).resp(id) + OTTF.popdata{j}.lfp.tfgamma;
                    avg(av).rate(oid) = oldrate;
                    avg(av).rate(id) = avg(av).rate(id) + OTTF.popdata{j}.lfp.rates./max(OTTF.popdata{j}.lfp.rates);
                    avg(av).alln(id) = avg(av).alln(id) + 1;
                end
                plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
                hold on;
            end
            OTTF.plot.labels(j) = 0;
        end
        if subtype ==1
            xlabel('Depth');
            ylabel('Coherency');
        elseif subtype == 2
            xlabel('Nspikes');
            ylabel('Coherency');
        elseif subtype == 3
            xlabel('Coherent F');
            ylabel('Coherency');
        elseif subtype == 4
            xlabel('Depth');
            ylabel('Stimulus Following');
        elseif subtype == 5
            xlabel('Coherency');
            ylabel('Stimulus following');
        end
        runlist('store',OTTF);
        GetFigure(OTTF.tag.extra);
        nid = find(avg(1).alln > mean(avg(1).alln)/2);
        hold off;
        plot(avg(1).alltfs,avg(1).resp./avg(1).alln);
        hold on;
        plot(avg(1).alltfs(nid),avg(1).rate(nid)./avg(1).alln(nid),'ro-');
        nid = find(avg(2).alln > mean(avg(2).alln)/2);
        plot(avg(2).alltfs(nid),avg(2).rate(nid)./avg(2).alln(nid),'go-');
    elseif isfield(OTTF,'xfunc')
        OTTF = feval(OTTF.xfunc,'DATA',OTTF,'popplot');
    end
    runlist('store',OTTF);
elseif strmatch(list,'nextmarked')
    it = findobj(gcf, 'Tag',OTTF.listtag);
    n = get(it, 'value');
    idx = find(OTTF.data.marked == 1);
    id = min(find(idx > n));
    if ~isempty(id)
        set(it, 'value',idx(id));
        OTTF = doentry(OTTF.listids(idx(id)),OTTF,'setid');
    end
elseif strmatch(list,'nexti')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n < length(OTTF.listids);
       n = n+1;
       set(it, 'value',n);
       OTTF = doentry(OTTF.listids(n),OTTF,'setid');
     end
elseif strmatch(list,'prev')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n > 1
       n = n-1;
       set(it, 'value',n);
       OTTF = doentry(OTTF.listids(n),OTTF,'setid');
       OTTF.id = n;
     end
elseif strmatch(list,'save')
    if ~isempty(varargin) & ischar(varargin{1})
        OTTF.statefile = varargin{1};
    else
    OTTF.statefile = SetStateName(OTTF);
    end
    if exist(OTTF.statefile,'file') %confirm overwrite
        resp = questdlg(sprintf('%s exists. Overwrite?',OTTF.statefile),'Confrim Overwrite','Yes','No','Yes');
        if strcmp(resp,'No')
            return;
        end
    end
    fprintf('Saving Data to %s\n',OTTF.statefile);
    BackupFile(OTTF.statefile,'print');
    save(OTTF.statefile,'OTTF');
elseif strmatch(list,'loadstate')
    if length(varargin) > 0 & ischar(varargin{1})
        OTTF.statefile = varargin{1};
    else
    OTTF.statefile = SetStateName(OTTF);
    end
    OTTF = LoadState(OTTF);
       runlist('store',OTTF);
elseif strmatch(list,'comment')
     OTTF.comments{OTTF.id} = get(findobj('Tag','FileComment'),'String');
     runlist('store',OTTF);
elseif strmatch(list,'update')
     top = num2str(OTTF.toplevel);
     OTTF.prefix = get(findobj('Tag','Prefix'),'String');
     OTTF.args = get(findobj('Tag','Args'),'String');
     OTTF.plot.legendpos = get(findobj('Tag','LegendPos'),'value') -1;
     OTTF.state.verbose = get(findobj('Tag','Verbose'),'value');
     OTTF.plot.showsdf = get(findobj(OTTF.toplevel,'Tag','Sdf'),'value');
     if ~isempty(findobj('Tag',OTTF.tag.select))
         OTTF.state.refit = get(findobj('Tag','Refit'),'value');
     if isempty(OTTF.state.refit)
         OTTF.state.refit = 0;
     end
     OTTF.state.rebuild = get(findobj('Tag',['Rebuild' top]),'value');
     if isempty(OTTF.state.rebuild)
         OTTF.state.rebuild = 0;
     end
     OTTF.plot.plotall = get(findobj('Tag','PlotAll'),'value');
     if isempty(OTTF.plot.plotall)
         OTTF.plot.plotall = 0;
     end
     OTTF.plot.forcezero = get(findobj('Tag','ForceZero'),'value');
     OTTF.plot.autoplot = get(findobj('Tag','AutoPlot'),'value');
     OTTF.plot.smoothing = get(findobj('Tag','Smooth'),'value');
     OTTF.plot.spikeadccorr = get(findobj('Tag','ADCcorr'),'value');
     OTTF.plot.triglfp = get(findobj('Tag','TrigLFP'),'value');
     OTTF.plot.plotem = get(findobj('Tag','PlotEM'),'value');
     OTTF.plot.nresample = getstr2num('Nresample');
     OTTF.plot.pval = getstr2num('Pcrit');
     OTTF.plot.minpsychtrials = getstr2num('Npsych');
     OTTF.plot.setmonkey = get(findobj('Tag','SetMonkey'),'value')-1;
     end

     runlist('store',OTTF);
     if isfield(OTTF,'xfunc')
        OTTF = feval(OTTF.xfunc,'DATA',OTTF,'update');
     elseif OTTF.plot.autoplot
         PlotCurrent;
     end
elseif strmatch(list,'Replay')
    id = OTTF.id;
    load(OTTF.fstrings{id});
    ReplaySpikes(Expt.Header.LstName);
    figure(findobj('Tag',OTTF.tag.top));
elseif strmatch(list,'setplot')
     SetPlotTypes;
     if OTTF.plot.autoplot
         PlotCurrent;
     end
elseif strmatch(list,'setentry')
    PlotCurrent;
elseif strmatch(list,'relist')
     listtype= get(findobj('Tag','ListType'),'value');
     relist(listtype,OTTF);
end


function UpdateWindows(OTTF)

it = findobj('Tag',OTTF.tag.select);
cm = findobj('Tag','FileComment','Parent',it);
if ~isempty(cm)
    set(cm,'String',OTTF.comments{OTTF.id});
end


function OTTF = LoadState(OTTF)

    toplevel =OTTF.toplevel;
    oldstate = OTTF.state; 
    oldplot = OTTF.plot;
    oldgui = OTTF.gui;
    menus = OTTF.menus;
    oldDATA = OTTF;
    
    tic; load(OTTF.statefile); toc
    f = fields(oldplot);
    for j = 1:length(f)
        if ~isfield(OTTF.plot,f{j})
            OTTF.plot.(f{j}) = oldplot.(f{j});
        end
    end
    f = fields(oldstate);
    for j = 1:length(f)
        if ~isfield(OTTF.state,f{j})
            OTTF.state.(f{j}) = oldstate.(f{j});
        end
    end
    if ~isfield(OTTF.plot,'condense')
        OTTF.plot.condense = 0;
    end
    if ~isfield(OTTF.plot,'flip')
        OTTF.plot.flip = 0;
    end
    OTTF.gui = oldgui;
    OTTF.toplevel = toplevel;
    OTTF.menus = menus;
    f = fields(oldDATA);
    for j = 1:length(f)
        if ~isfield(OTTF,f{j})
            OTTF.(f{j}) = oldDATA.(f{j});
        end
    end
    if ~isempty(OTTF.xfunc)
        OTTF = feval(OTTF.xfunc, 'checkload', OTTF);
    end

function name = SetStateName(OTTF)

if isempty(OTTF.statefile)
    if strfind(OTTF.list,'.mall')
        OTTF.statefile = strrep(OTTF.list,'.mall','.mat');
    else
        OTTF.statefile = [OTTF.list '.mat'];
    end
end
name = OTTF.statefile;

function RePlot(OTTF)

if isfield(OTTF,'xfunc')
    feval(OTTF.xfunc,'popplot');
else
    GetFigure(OTTF.tag.extra);
    hold off;
    for j = 1:length(OTTF.plot.allpts)
        plot(OTTF.plot.allpts(1,j),OTTF.plot.allpts(2,j),'o','buttondownfcn',['runlist(''TouchPoint'',gcf,' num2str(j) ')']);
        hold on;
    end
end

function OTTF = doentry(id, OTTF,varargin)
   
j = 1;
while(j < nargin - 1)
    if strncmpi(varargin{j},'setid',5)
        OTTF.id = id;
        set(OTTF.toplevel,'Name',sprintf('%s %d/%d',OTTF.list,id,length(OTTF.fstrings)));
    end
    j = j+1;
end

funcalled = 0;
argon = {};
  if get(findobj(OTTF.toplevel,'Tag','LFP'),'value')
       argon = {argon{:} {'lfp'}};
   end
   if get(findobj(OTTF.toplevel,'Tag','ShowN'),'value')
       argon = {argon{:} {'ShowN'}};
   end
   if get(findobj(OTTF.toplevel,'Tag','Pcolor'),'value')
       argon = {argon{:} {'pcolor'}};
   end
   if OTTF.plot.plotall
       argon = {argon{:} {'plotall'}};
   end
   if OTTF.plot.showsdf
       argon = {argon{:} {'sdfall'}};
   end
   if length(OTTF.args) > 1
       commas = strfind(OTTF.args,',');
       if isempty(commas)
           argon = {argon{:} OTTF.args};
       else
           last = 1;
           for j = 1:length(commas)
               argon = {argon{:} OTTF.args(last:commas(j)-1)};
               last = commas(j)+1;
           end
           argon = {argon{:} OTTF.args(commas(j)+1:end)};
       end
   end
   GetFigure(OTTF.tag.fig);
 % if xfunc is defined, it must handle the Spike replay
 if strmatch(OTTF.plot.onetype,'Spike') & ~isfield(OTTF,'xfunc')
     load(OTTF.fstrings{id});
     if exist('cExpt','var')
         PlotExptSpikes(cExpt);
%         if isfield(DATA,'dprime')
%             OTTF.popdata.spkdprime(id) = DATA.dprime;
%         end
     else
     CheckSpike(Expt);
     end
 elseif strmatch(OTTF.plot.onetype,'Replay')
     if exist(OTTF.fstrings{id},'file')
         res = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,argon{:});
         if ~isempty(res.Data)
             ReplaySpikes(res.Data.Header.LstName);
             figure(findobj('Tag',OTTF.tag.top));
         else
             fprintf('Error: No expt in %s\n',OTTF.fstrings{id});
         end
     end
 elseif strmatch('FitSine',OTTF.plot.onetype)
     subplot(2,1,1);
     hold off;
     data = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,'noline',argon{:});
     hold on;
     for j = 1:length(data.linevals);
         fits(j) = FitSine(data.x(:,j),data.means(:,j));
         fits(j).fitx = [min(data.x(:,j)):0.1:max(data.x(:,j))];
         fits(j).fitted = FitSine(fits(j).fitx,fits(j).answer,'eval');
         plot(fits(j).fitx,fits(j).fitted,'color',data.colors{j});
     end
     subplot(2,1,2);
     hold off;
     params = [fits.answer];
     plot(data.linevals,params(1,:),'o');
     hold on;
     plot(data.linevals,abs(params(2,:)),'ro');
 elseif strmatch(OTTF.plot.onetype,'Psych')
     data = PlotExpt(OTTF.fstrings{id},argon{:},'psych'); 
     if isfield(data,'conc')
         OTTF.popdata{id}.lfp.conc = res.conc(2);
     end
 elseif strmatch(OTTF.plot.onetype,'Consec')
     data = PlotExpt(OTTF.fstrings{id},'consec','Pd',argon{:});
     for j = 1:length(data)
         ids = find(data(j).n(:) > 0);
         ncount(j) = sum(data(j).n(ids));
         sumcount(j) = sum(data(j).means(ids));
     end
     if length(data) > 3 & ~isempty(data(2).psych) & ~isempty(data(4).psych)
         OTTF.popdata{id}.psf(1) = data(2).psych.fit(2);
         OTTF.popdata{id}.psf(2) = data(4).psych.fit(2);
         OTTF.popdata{id}.ntrials = sum([data(4).psych.data.n]) + sum([data(2).psych.data.n]);
         OTTF.popdata{id}.mean(1) = sum(sumcount(1:2))./sum(ncount(1:2));
         OTTF.popdata{id}.mean(2) = sum(sumcount(3:4))./sum(ncount(3:4));
         fprintf('%s dTh = %.4f\n',data(1).name,diff(OTTF.popdata{id}.psf));
     end
 elseif strmatch(OTTF.plot.onetype,'VHDisp')
     data = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,argon{:});
     theta = pi/2 - (GetEval(data.Data,'Ro') * pi/180);
     xi = data.x;
     yi = data.y;
     if OTTF.plot.smoothing
         xi = linspace(min(data.x),max(data.x),range(data.x)/40);
         yi = linspace(min(data.y),max(data.y),range(data.y)/40);
         [xxi, yyi] = meshgrid(xi,yi);
         xinc = 0.1;
         yinc = 0.1;
         %Gauassian smoothing interpolation instead
         if(xinc < yinc)
             smoothing = yinc/2;
         else
             smoothing = xinc/2;
         end
         zi = Interpf(X,Y,Z,xxi,yyi,1,smoothing);
     end
         data.dx = data.x .* cos(theta) + data.y * sin(theta);
         data.dy = data.y * cos(theta) - data.x * sin(theta);
         [X,Y,Z] = fillpmesh(data.dx,data.dy,data.means);
         hold off;
     pcolor(X,Y,Z);
     if OTTF.plot.forcezero
         zlim = caxis;
         caxis([0 zlim(2)]);
     end
     colorbar;
     axis('image');
     StoreData(OTTF,data,id);
     title(splitpath(data.name));
 elseif strmatch(OTTF.plot.onetype,'Spatial')
     GetFigure(OTTF.tag.fig);
     datafile = strrep(OTTF.fstrings{id},'rds','grating');
     [s, t] = regexp(datafile,'\.[A-z]*[A-Z]\.');
     if isempty(s)
         return;
     end
     suffix = datafile(s(1)+1:t(1)-1);
     subplot(2,2,1);
     tfname = strrep(datafile,suffix,'TF');
     if ~exist(tfname,'file')
         tfname = strrep(datafile,suffix,'TFM');
     end
     delete(allchild(gca));
     res = PlotExpt(tfname);
     if isempty(res)
         title(['No ' splitpath(tfname)]);
     end
     subplot(2,2,2);
     sfname = strrep(datafile,suffix,'SF');
     delete(allchild(gca));

     if ~exist(sfname,'file')
         sfname = strrep(datafile,suffix,'SFM');
     end

     if isempty(strfind(sfname,suffix))
         res = PlotExpt(sfname);
         if isempty(res)
             title(['No ' splitpath(sfname)]);
         end
     end
     subplot(2,2,3);
     otname = strrep(datafile,suffix,'OXM');
     delete(allchild(gca));

     if ~exist(otname,'file')
         otname = strrep(datafile,suffix,'OT');
     end

     if isempty(strfind(otname,suffix))
         res = PlotExpt(otname);
         if isempty(res)
             title(['No ' splitpath(otname)]);
         end
     end
     if isfield(OTTF,'xfunc')
         OTTF = feval(OTTF.xfunc,'DATA',OTTF,'Tuning');
     end
 elseif strmatch(OTTF.plot.onetype,'SpikeTwo','exact')
     GetFigure(OTTF.tag.fig);
     load(OTTF.fstrings{id});
     bExpt = Expt;
     clfirst = strrep(OTTF.fstrings{id},'c2','c1');
     load(clfirst);
     CheckSpikePair(Expt,bExpt,'dispw',1000);
 elseif strmatch(OTTF.plot.onetype,'SpikeW','exact')
     GetFigure(OTTF.tag.fig);
     subplot(1,1,1);
     load(OTTF.fstrings{id});
     ds = SpikeShape(Expt.Spike,'plot');
 elseif strmatch(OTTF.plot.onetype,'TwoCluster')
          subplot(1,1,1);
         thedata = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,argon{:});
         clfirst = strrep(OTTF.fstrings{id},'c2','c1');
         thedata = PlotExpt(clfirst,'legendpos',OTTF.plot.legendpos,argon{:});
         thedata = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,argon{:},'hold','Nlin',1);
 elseif strmatch(OTTF.plot.onetype,{'Collapse1'})
       thedata = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,'Collapse',1,argon{:});
    
 elseif strmatch(OTTF.plot.onetype,{'Tuning' 'TrigLFP'})
          subplot(1,1,1);
     if(OTTF.state.verbose)
         tic;
     end
     if strmatch(OTTF.plot.onetype,'TrigLFP')
         argon = {argon{:} 'lfptrig'};
     end
     if OTTF.plot.condense
         argon = {argon{:} 'condense'};
     end
     if OTTF.state.usedefaultplot
         thedata = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,argon{:});
         plotdata = StoreData(OTTF,thedata,id);
         OTTF.data.current = thedata;
         if ~isfield(thedata,'title') | isempty(thedata(1).title)
             if isfield(thedata,'name')
             title(splitpath(thedata.name));
             end
         end
         if isfield(thedata,'name') & strfind(thedata.name,'sqcorrug.OT');
             OTTF.dm(id) = GetEval(thedata.Data,'dm');
             OTTF.meandisp(id) = GetEval(thedata.Data,'dx');
         end
     end
     if(OTTF.state.verbose)
         toc
     end
     if OTTF.state.condense
         data = Condense(plotdata{id});
         fprintf('%.2f,%.2f  %.2fx%.2f\n',data.x,data.y,data.width,data.height);
     end
     if isfield(OTTF,'xfunc')
         OTTF = feval(OTTF.xfunc,'DATA',OTTF,'Tuning');
         funcalled = 1;
     end
 elseif isfield(OTTF,'xfunc')

    OTTF = feval(OTTF.xfunc,'DATA',OTTF,varargin{:});
     funcalled = 1;
     OTTF.test = 1;
 else
     subplot(1,1,1);
     data = PlotExpt(OTTF.fstrings{id},'legendpos',OTTF.plot.legendpos,argon{:});
     StoreData(OTTF,data,id);
     if(OTTF.state.verbose)
         data = Condense(OTTF.data.plotdata{id});
         fprintf('%.2f,%.2f  %.2fx%.2f\n',data.x,data.y,data.width,data.height);
     end
 end

 if strmatch(OTTF.plot.extraplot,'SpikeTwo')
     GetFigure(OTTF.tag.extra);
     argon = {};
     load(OTTF.fstrings{id});
     bExpt = Expt;
     if OTTF.plot.spikeorr
         argon = {argon{:}, 'adcc'};
     end
     clfirst = strrep(OTTF.fstrings{id},'c2','c1');
     load(clfirst);
     CheckSpikePair(Expt,bExpt,'dispw',1000,argon{:});
 elseif strmatch(OTTF.plot.extraplot,'Tuning')
     GetFigure(OTTF.tag.extra);
     argon = {};
     load(OTTF.fstrings{id});
     PlotExpt(Expt);
 elseif strmatch(OTTF.plot.extraplot, OTTF.plot.addextra) & isfield(OTTF,'xfunc')
 end

 subtype = strmatch(OTTF.plot.xop,{'TFLfp','SZLfp','TFLFPall'});
 if  subtype & ~isempty(thedata) & isfield(thedata,'lfpwr')
     GetFigure(OTTF.tag.fig);
     thedata.Data.Header.lfpfrq = thedata.lfpfrq;
     if subtype == 2
         res = TFFuncs(OTTF, thedata.Data, thedata, 'bysize');
     elseif subtype == 3
         res = TFFuncs(OTTF, thedata.Data, thedata, 'tfcoh');
     else
         res = TFFuncs(OTTF, thedata.Data, thedata);
     end
     OTTF.popdata{id}.lfp.conc = res.conc(2);
     OTTF.popdata{id}.lfp.coh_peakf = res.peakf;
     OTTF.popdata{id}.lfp.nspk = res.nspk;
     OTTF.popdata{id}.lfp.tfs = res.tfs;
     OTTF.popdata{id}.lfp.tfgamma = res.clfp;
     OTTF.popdata{id}.lfp.rates = res.rates;
     OTTF.popdata{id}.lfp.stimcoh = res.fpwrfrac;
     [area, depth] = GetPenData(OTTF.data.current.name);
     OTTF.popdata{id}.area = area;
     if isempty(depth)
         OTTF.popdata{id}.depth = -4000;
     else
         OTTF.popdata{id}.depth = depth;
     end
 elseif strmatch(OTTF.plot.xop, OTTF.plot.addxop) & isfield(OTTF,'xfunc')
  %      OTTF = feval(OTTF.xfunc,'DATA',OTTF,'xop',varargin{:});
 end
 
if isfield(OTTF,'xfunc') & ~funcalled
% it is up to xfunc to make sure that all three plot types are handled, and
% nothin is done if nothing is necessary
     OTTF = feval(OTTF.xfunc,'DATA',OTTF,varargin{:});
     
 end
runlist('store',OTTF);
if OTTF.state.verbose & isfield(OTTF.data, 'current')
    [area, depth] = GetPenData(OTTF.data.current.name);
    fprintf('%s V%d %.3f\n',splitpath(OTTF.data.current.name),area,depth);
    toc;
end


UpdateWindows(OTTF);


function PlotCurrent()
 [OTTF, n] = runlist('getstate');

  doentry(OTTF.listids(n),OTTF,'setid');

function SetPlotTypes()
%SetPlotTypes looks at the values of the menus for selecting
%plot types, and records them in OTTF for use elsewhere

  OTTF = runlist('getstate');

  it = findobj(OTTF.toplevel,'Tag','plottype');
  type = get(it, 'value');
  str = get(it, 'String');
  OTTF.plot.onetype = strrep(str(type,:),' ',''); %The string shown on the menu.
 
  it = findobj('Tag','extraplot');
  type = get(it, 'value');
  str = get(it, 'String');
  OTTF.plot.extraplot = strrep(str(type,:),' ',''); %The string shown on the menu.
 
  it = findobj('Tag','xop');
  if ~isempty(it)
      type = get(it, 'value');
      str = get(it, 'String');
      OTTF.plot.xop = strrep(str(type,:),' ',''); %The string shown on the menu.
  end
  
  it = findobj('Tag','sizetype');
  type = get(it, 'value');
  OTTF.plot.sizetype = type; 
  runlist('store',OTTF);
  doentry(OTTF.id,OTTF);


function LoadOTTFData(OTTF, j)
if isfield(OTTF,'xfunc')
    OTTF = feval(OTTF.xfunc,'LoadOne');
    return;
end
        
function OTTF = LoadAllData(OTTF)

    if isfield(OTTF,'xfunc')
        OTTF = feval(OTTF.xfunc,'LoadAll');
        return;
    end

for j = 1:length(OTTF.fstrings)
  fprintf('%s',OTTF.fstrings{j});
  LoadOTTFData(OTTF, OTTF.fstrings{j});
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


function OTTF = relist(step,OTTF)

n = 1;
strs = get(findobj('Tag','ListType'),'String');
OTTF.listtype = strs{get(findobj('Tag','ListType'),'value')};
OTTF.listids = [];

if strmatch(OTTF.listtype,'All')
    lst = findobj('Tag',OTTF.listtag);
  set(lst, 'String',OTTF.fstrings);
  OTTF.listids = [1:length(OTTF.fstrings)];
  runlist('store',OTTF);
  return;
end
if isfield(OTTF,'xfunc')
%Call xfunc to set list order.
    OTTF = feval(OTTF.xfunc,'DATA',OTTF,'setlistorder');
    if isfield(OTTF,'listorder')
        OTTF.listids = OTTF.listorder;
        for j = 1:length(OTTF.listids)
            fstrings{j} = OTTF.fstrings{OTTF.listids(j)};
        end
    end
end

if isempty(OTTF.listids)
    if strmatch(OTTF.listtype,{'rufus', 'dufus', 'lem', 'ic', 'dae', 'jbe'})
        for j = 1:length(OTTF.fstrings)
            if ~isempty(strfind(OTTF.fstrings{j},OTTF.listtype))
                fstrings{n} = OTTF.fstrings{j};
                OTTF.listids(n) = j;
                n = n+1;
            end
        end
    end
else
    for j = 1:length(OTTF.fstrings)
        load(OTTF.fstrings{j});
        sf = GetEval(Expt,'sf','mode');
        fs = GetEval(Expt,'Fs','mode');
        OTTF.sfsteps(j) = fs * sf;
        x = mod(OTTF.sfsteps(j),1.0);
        if((x < 0.9 & x > 0.1 & step == 0.5) | step == 0 | ( x < 0.1 & step ==1))
            fstrings{n} = OTTF.fstrings{j};
            OTTF.listids(n) = j;
            n = n+1;
        end
    end
end
  lst = findobj('Tag',OTTF.listtag);
  set(lst, 'String',fstrings);
  runlist('store',OTTF);
% InitOTTF makes the basic gui
%
%
%
%
  
  
function OTTF = InitOTTF(list, OTTF, varargin)
global bgcfileprefix;

nogui = 0;
wsize = [27 35]; %in columns
if iscellstr(list)
    OTTF.fstrings = list;
elseif ~exist(list,'file')
    fprintf('No File %s',list);
    return;
else
    OTTF.fstrings = textread(list,'%s');
end

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
OTTF.plot.showsdf = 0;
OTTF.state.storedata = 0;
OTTF.state.condense = 0;
OTTF.state.includeratio = 0.5;
OTTF.state.usedefaultplot = 1;
OTTF.state.dirw = 0.5;
OTTF.state.refit = 0;

OTTF.plot.plotall = 0;
OTTF.plot.onetype = 'Tuning';
OTTF.plot.xop = 'none';
OTTF.plot.addxop = {};
OTTF.plot.addextra = {};
OTTF.plot.sizetype = 1;
OTTF.plot.poptype = 1;
OTTF.plot.sdfsmooth = 150;
OTTF.plot.sdftype = 'exp';
OTTF.plot.legendtype = 1;
OTTF.plot.sdfstart = -5000;
OTTF.plot.sdfend = 5000;
OTTF.plot.sdfstep = 20;
OTTF.plot.legendpos = 2;
OTTF.plot.extraplot = 'None';
OTTF.plot.forcezero = 0;
OTTF.plot.autoplot = 0;
OTTF.plot.smoothing = 0;
OTTF.plot.spikeadccorr = 0;
OTTF.plot.triglfp = 0;
OTTF.plot.plotem = 0;
OTTF.plot.nresample = 0;
OTTF.plot.minpsychtrials = 0;
OTTF.plot.pval = 0.05;
OTTF.plot.setmonkey = 0;
OTTF.plot.plotem = 0;
OTTF.plot.condense = 0;
OTTF.plot.setmonkey = 0;
OTTF.state.nogui = 0;
OTTF.statefile = [];
OTTF.args = '';


j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'sizestrings',7)
        OTTF.sizestrings = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'size',2)
        j = j+1;
        wsize = varargin{j};
    elseif strncmpi(varargin{j},'nogui',5)
        j = j+1;
        nogui = 1;
    end
    j = j+1;
end

scrsz = get(0,'Screensize');

if scrsz(3) == 1920 & scrsz(4) == 1200
    wsc = 1.1;
elseif scrsz(3) == 1024 & scrsz(4) == 768
    wsc = 0.9;
else
    wsc = 1;
end

OTTF.wsc = wsc;
cw = 10 * wsc;
bh = 9 * wsc;
listlen = 15;
wsiz = [wsize(1)*bh cw * wsize(2)];
for j = 1:length(OTTF.fstrings)
    if OTTF.fstrings{j}(1) == '#';
%        delete(OTTF.fstrings{j});
    end
end
OTTF.gui.ch = bh;
OTTF.gui.cw = cw;
OTTF.gui.wsiz = wsiz;

if ~isempty(bgcfileprefix) & ~strncmp(bgcfileprefix,OTTF.fstrings{1},length(bgcfileprefix))
    OTTF.prefix = bgcfileprefix;
else
    OTTF.prefix = '';  
end
%if use /bgc/bgc as path, should always work.
if strfind(OTTF.fstrings{1},'/bgc/bgc')
OTTF.prefix = '';
else
OTTF.prefix = '';
end
if iscellstr(list)
    OTTF.list = 'struct';
else
OTTF.list = list;  
end

j = 1;
while(j < nargin-1)
    if strncmpi(varargin{j},'prefix',5)
        j = j+1
        OTTF.prefix = varargin{j};
    elseif strncmpi(varargin{j},'verbose',5)
        j = j+1
        OTTF.state.verbose = 1;
    elseif strncmpi(varargin{j},'wsiz',4)
        j = j+1;
        wsiz = varargin{j};
    elseif strncmpi(varargin{j},'state',5)
        j = j+1;
        OTTF.statefile = varargin{j};
        if strncmp(varargin{j-1},'State',5)
            OLDDAT = OTTF;
            tic;
            load(varargin{j});
            toc;
            OTTF.tag = OLDDAT.tag;
            OTTF.statefile = OLDDAT.statefile;
        end
    
    end
    j = j+1;
end


if ~iscellstr(list)
OTTF.fstrings = ReadCellList(list);
end
OTTF.listids = 1:length(OTTF.fstrings);
laststr = 0;
for j = 1:length(OTTF.fstrings)
    OTTF.fstrings{j} = strrep(OTTF.fstrings{j},'\','/');
    if OTTF.fstrings{j}(1) == '#';
        deletestr(j) = 1;
        laststr = 1;
    elseif laststr & isempty(strfind(OTTF.fstrings{j},'.mat')) %% still in comment
        deletestr(j) = 1;
    else
        goodstr(j) = 1;
        laststr = 0;
    end
end
OTTF.fstrings = sort({OTTF.fstrings{find(goodstr)}});


for j = 1:length(OTTF.fstrings)
    if ~isfield(OTTF,'comments') | j > length(OTTF.comments) | isempty(OTTF.comments{j})
        OTTF.comments{j} = splitpath(OTTF.fstrings{j});
    end
end

if ~isempty(OTTF.prefix)
    for j = 1:length(OTTF.fstrings)
        OTTF.strings{j} = OTTF.fstrings{j};
        OTTF.fstrings{j} = [OTTF.prefix OTTF.fstrings{j}];
    end
end

for j = 1:length(OTTF.fstrings)
    OTTF.data.marked(j) = 0;
end


if ischar(list)
    idx = findstr(list,'/');
if isempty(idx)
  idx = findstr(list,'\\');
end
if isempty(idx)
  listname = list;
else
  listname = list(idx(end)+1:end);
end
else
    listname = 'struct';
end

if nogui
    OTTF.state.nogui = 1;
    return;
end

  cntrl_box = figure('Position', [100 scrsz(4)-(wsiz(2)+bh*4) wsiz(1) wsiz(2)],...
       'NumberTitle', 'off', 'Tag',OTTF.tag.top,'Name',sprintf('OTTF %s (%d)',listname,length(OTTF.fstrings)));
  lst = uicontrol(gcf, 'Style','listbox','String',OTTF.fstrings,...
		'Callback', ' runlist(''setentry'')','Tag',OTTF.listtag,...
		'Position',[10 10 wsiz(1)-20 cw*listlen]);

  SPACE = 2 * wsc;
  VSPACE = 5 * wsc;

  OTTF.toplevel = gcf;
  OTTF.gui.lst = lst;
  bp(1) = SPACE; bp(2) = cw*(listlen+1); bp(3) = 40; bp(4) = 22;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runlist(''nexti'')',...
'String', '>>', 'Tag',OTTF.tag.data,'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runlist(''prev'')',...
'String', '<<', 'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runlist(''All'')',...
'String', 'All', 'Position', bp);

bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runlist(''More'')',...
'String', 'More', 'Position', bp);


bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' runlist(''NewPrefix'')',...
'String', OTTF.prefix, 'Tag', 'Prefix','Position', bp);
  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' runlist(''Replay'')',...
'String', 'Replay', 'Position', bp,'Tag','Replay');


    bp(1) = SPACE; bp(3) = 25; bp(2) = bp(2) + bp(4)+VSPACE; bp(4) = 22;
bp(3) = 50;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'LFP', 'Tag', 'LFP', 'Position', bp);

  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw*6;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'ShowN', 'Tag', 'ShowN', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 6*cw;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Pcolor', 'Tag', 'Pcolor', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runlist(''update'')',...
'String', 'sdf', 'Tag', 'Sdf', 'Position', bp);
bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' runlist(''update'')',...
'String', 'Verbose', 'Tag', 'Verbose', 'Position', bp);

  bp(1) = SPACE;
  bp(2) = bp(2) +bp(4) + VSPACE;
  bp(3) = 70;
    uicontrol(gcf,'style','pop','string','Tuning|Collapse 1|Collapse 2|Spike|Spike Two|SpikeW|Replay|Two Cluster|TrigLFP|VHDisp|OTuning|FitSine|Spatial|Psych|Consec', ...
		    'Callback', ' runlist(''setplot'')', 'Tag','plottype',...
		    'position',bp,'value',1);

   bp(1) = bp(1)+bp(3)+SPACE;
   bp(3) = 30;
  uicontrol(gcf,'Style', 'text','String','Extra','Position', bp);        
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 50;
  uicontrol(gcf,'style','pop','string','None|Tuning|Spike|Spike Two|sacortho|Jump|SacJump|spike|TF|SZ|OP|OT|SF|SZcmp|Psych', ...
		    'Callback', ' runlist(''setplot'')', 'Tag','extraplot',...
		    'position',bp,'value',1);
 
   bp(1) = bp(1)+bp(3)+SPACE;
   bp(3) = 20;
  uicontrol(gcf,'Style', 'text','String','Op','Position', bp);   
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 50;
  uicontrol(gcf,'style','pop','string','None|TFLfp|TFOpt|TFgammaLFP|SZLfp|Combine|PJump|OJump', ...
		    'Callback', ' runlist(''setplot'')', 'Tag','xop',...
		    'position',bp);
  bp(1) = SPACE;
  bp(2) = bp(2) + bp(4) + VSPACE;
  bp(3) = 40;
  if isfield(OTTF,'sizestrings')
  uicontrol(gcf,'Style', 'text','String','Scheme','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
      uicontrol(gcf,'style','pop','string',OTTF.sizestrings, ...
		    'Callback', ' runlist(''setplot'')', 'Tag','sizetype',...
		    'position',bp);
  else
  uicontrol(gcf,'Style', 'text','String','Size','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
      uicontrol(gcf,'style','pop','string','All|First|Last|Both|Combine|2|3|4|5|6', ...
		    'Callback', ' runlist(''setplot'')', 'Tag','sizetype',...
		    'position',bp);
  end
  bp(1) = bp(1) +bp(3) + SPACE;
  uicontrol(gcf,'Style', 'text','String','Legend','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','Auto|UR|UL|LL|LR|Off|None', ...
		    'Callback', ' runlist(''update'')', 'Tag','LegendPos',...
		    'position',bp);
  
  bp(1) = SPACE;
  bp(2) = bp(2) + bp(4) + VSPACE;
  bp(3) = cw * 4;
  uicontrol(gcf,'Style', 'text','String','List','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','All|FullPeriod|HalfPeriod|dufus|rufus|lem|ic|dae|jbe', ...
		    'Callback', ' runlist(''relist'')', 'Tag','ListType',...
		    'position',bp);
  
  bp(1) = bp(1)+bp(3) +SPACE;
  uicontrol(gcf,'Style','pushbutton', 'Callback', ' runlist(''popplot'')','String','PopPlot','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','Isolation|LFPCoh|LFPCohF|LFPnspk|LFPstim|Psych||ConsecMeans|ConsecLFP|SpikeW|New',...
		    'Callback', ' runlist(''popplot'')', 'Tag','PopPlot',...
		    'position',bp);

  bp(1) = SPACE;
  bp(2) = bp(2) + bp(4) + VSPACE;
  bp(3) = cw*4+SPACE;
  uicontrol(gcf,'Style', 'text','String','Args','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*40;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' runlist(''update'')',...
      'String', OTTF.args, 'Tag', 'Args','Position', bp);
        
  x = 10; y = 180;

  hm = uimenu(gcf,'Label','Mark');
  uimenu(hm,'Label','&Mark','Callback','runlist(''Mark'')');
  uimenu(hm,'Label','UnMark','Callback','runlist(''UnMark'')');
  uimenu(hm,'Label','List','Callback','runlist(''PrintMark'')');
  uimenu(hm,'Label','UnMark All','Callback','runlist(''UnMarkAll'')');
  uimenu(hm,'Label','Load One','Callback',' runlist(''loadone'')');
  uimenu(hm,'Label','Load All','Callback',' runlist(''load'')');
  uimenu(hm,'Label','Save State','Callback',' runlist(''save'')');
  uimenu(hm,'Label','Save Selected','Callback',{@SaveList, 'selected'});
  uimenu(hm,'Label','Load State','Callback',' runlist(''loadstate'')');
  uimenu(hm,'Label','Close','Callback',' runlist(''close'')');
  OTTF.menus.mark = hm;
  hm = uimenu(gcf,'Label','Options','Tag','Optionmenu');
  uimenu(hm,'Label','Options','Callback',' runlist(''popselect'')');
  uimenu(hm,'Label','Strip list prefix','Callback',' runlist(''stripprefix'')');
  OTTF.menus.options = hm;
  set(gcf,'Menubar','none');
  hm = uimenu(gcf,'Label','&Next','Callback','runlist(''nextmarked'');');
  OTTF.menus.next = hm;

OTTF.figid.mainplot = figure('Tag',OTTF.tag.fig,'Name',OTTF.tag.fig,'Position',[300 scrsz(4)-470 512 ...
		    400]);
OTTF.id = 1;
OTTF.statefile = SetStateName(OTTF);
set(cntrl_box,'UserData',OTTF);


function SaveList(a, b, fcn)

DATA = GetDataFromFig(a);

if isfield(DATA.plot,'selected')
    id = find(DATA.plot.selected);
    str = DATA.fstrings(DATA.listids(id));
    name = uiputfile('./newlist.mall');
    if name
    fid = fopen(name,'w');
    for j = 1:length(str)
        fprintf(fid,'%s\n',str{j});
    end
    fclose(fid);
    end
end


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
runlist('store',OTTF);


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

 
  wsc = OTTF.wsc;
  SPACE = 3 * wsc;
  VSPACE = 5 * wsc;
  h = 220 * wsc;
  w = 350 * wsc;
  scrsz = get(0,'Screensize');
 
  cntrl_box = figure('Position', [200 scrsz(4)-(h+30)*wsc w*wsc h*wsc], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag',tag,'Name','Section Criteria');
   OTTF.figid.options = cntrl_box;
   top = num2str(OTTF.toplevel); 
   cw = 10 * wsc;
   ch = 11 * wsc;
   bh = 18*wsc;
   bp(1) = SPACE;
   bp(2) = ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
  
  uicontrol(gcf,'Style', 'CheckBox','String','Force Zero','Position', bp,...
      'Tag','ForceZero','Callback','runlist(''update'')','value',OTTF.plot.forcezero);
  bp(1) = bp(1) + bp(3) + SPACE;
  bp(3) = cw*7;
    uicontrol(gcf,'Style', 'CheckBox','String','Smooth','Position', bp,...
      'Tag','Smooth','Callback','runlist(''update'')','value',OTTF.plot.smoothing);
  bp(1) = bp(1) + bp(3) + SPACE;
  bp(3) = cw * 5;
  uicontrol(gcf,'Style', 'CheckBox','String','Auto','Position', bp,...
      'Tag','AutoPlot','Callback','runlist(''update'')','value',OTTF.plot.autoplot);
  bp(1) = bp(1) + bp(3) + SPACE;
  bp(3) = cw*7;
  uicontrol(gcf,'Style', 'CheckBox','String','Plot All','Position', bp,...
      'Tag','PlotAll','Callback','runlist(''update'')','value',OTTF.plot.plotall);
  bp(1) = bp(1) + bp(3) + SPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Refit','Position', bp,...
      'Tag','Refit','Callback','runlist(''update'')','value',OTTF.state.refit);
  bp(1) = bp(1) + bp(3) + SPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Rebuild','Position', bp,...
      'Tag',['Rebuild' top],'Callback','runlist(''update'')','value',OTTF.state.rebuild);
  bp(1) = SPACE;
  bp(2) = bp(2) + ch+ VSPACE;  
  bp(3) = w*wsc - bp(1);
  uicontrol(gcf,'Style', 'edit', 'Callback', ' runlist(''comment'');',...
    'String', OTTF.comments{OTTF.id}, 'Tag', 'FileComment','Position', bp);

  bp(1) = SPACE;
  bp(2) = bp(2) + ch+ VSPACE;  
  bp(3) = cw*7;
  uicontrol(gcf,'Style', 'text','String','Nresample','Position', bp);
  bp(1) = bp(1) + SPACE + bp(3);
  bp(3) = cw *5;
  uicontrol(gcf,'Style', 'Edit','String','Nresample','Position', bp,...
      'Tag','Nresample','Callback','runlist(''update'')','String',num2str(OTTF.plot.nresample));

  bp(1) = bp(1) + bp(3) + SPACE;
  uicontrol(gcf,'Style', 'text','String','Pcrit','Position', bp);
  bp(1) = bp(1) + SPACE + bp(3);
  bp(3) = cw *5;
  uicontrol(gcf,'Style', 'Edit','Position', bp,...
      'Tag','Pcrit','Callback','runlist(''update'')','String',num2str(OTTF.plot.pval));

  bp(1) = bp(1) + bp(3) + SPACE;
  uicontrol(gcf,'Style', 'text','String','Npsych','Position', bp);
  bp(1) = bp(1) + SPACE + bp(3);
  bp(3) = cw *5;
  if ~isfield(OTTF.plot,'minpsychtrials')
      OTTF.plot.minpsychtrials = 100;
  end
  uicontrol(gcf,'Style', 'Edit','Position', bp,...
      'Tag','Npsych','Callback','runlist(''update'')','String',num2str(OTTF.plot.minpsychtrials));

  
  if isempty(OTTF.plot.spikeadccorr)
      OTTF.plot.spikeadccorr = 0;
  end
  bp(1) = SPACE;
  bp(2) = bp(2) + ch*2;  
  bp(3) = cw*7;
  uicontrol(gcf,'Style', 'CheckBox','String','ADCcorr','Position', bp,...
      'Tag','ADCcorr','Callback','runlist(''update'')','value',OTTF.plot.spikeadccorr);
  bp(1) = bp(1)+bp(3) + SPACE;
  bp(3) = cw*6;
  if isempty(OTTF.plot.triglfp)
      OTTF.plot.triglfp = 0;
  end
      uicontrol(gcf,'Style', 'CheckBox','String','TrigLFP','Position', bp,...
      'Tag','TrigLFP','Callback','runlist(''update'')','value',OTTF.plot.triglfp);
  if ~isfield(OTTF.plot,'plotem')
      OTTF.plot.plotem = 0;
  end
  bp(1) = bp(1)+bp(3) + SPACE;
  bp(3) = cw*6;
      uicontrol(gcf,'Style', 'CheckBox','String','EyePos','Position', bp,...
      'Tag','PlotEM','Callback','runlist(''update'')','value',OTTF.plot.plotem); 
  bp(1) = bp(1)+bp(3) + SPACE;
  bp(3) = cw*6;
      uicontrol(gcf,'Style', 'pop','String','All|Lem|ruf|duf|ica','Position', bp,...
      'Tag','SetMonkey','Callback','runlist(''update'')','value',OTTF.plot.setmonkey+1); 
  OTTF.wsc = wsc;
  OTTF.cntrlbox = cntrl_box;
  OTTF.selectlast = bp(2);
  runlist('store',OTTF);
else
  cntrl_box = findobj('Tag',tag);
  figure(cntrl_box);
end
