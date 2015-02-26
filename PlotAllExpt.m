function PlotAllExpt(A, varargin)
%plots data from an AllExpt Struct
%takes a struct with all trials so far in a file, and plots most recent
%data. Intended to replace the plots in binoc eventually.

if isfield(Expt,'Expts') && isfield(Expt.Expts,'Start')
    expid = length([Expt.Expts.Start]);
else
    expid = [];
end
dattime = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'Expts',4)
        j = j+1;
        expid = varargin{j};
    elseif strncmpi(varargin{j},'NewData',4)
        dattime = now;
    end
    j = j+1;
end

[ofig, isnew] = GetFigure('OnlinePlot');
if isnew
    cw = 8;
    ch = 12;
    DATA.plot.probe = 1;
    DATA.plot.shown = 1;
    DATA.plot.psych = 0;
    DATA.toplevel = ofig;
    DATA.probe = 1;
    DATA.nprobes = 8;
    DATA.cluster = 2;
    DATA.dattime = 0;
    DATA.ncalls = 0;
    bp(1) = 10; bp(2) = 10;
    bp(3) = cw * 4; bp(4) = ch;
    uicontrol(ofig,'style','pop','string',sprintf('%d',[1:DATA.nprobes])','Position',bp,'Tag','SetProbe',...
        'Callback',@Update);
     bp(1) = bp(1) + bp(3) + 10;
    uicontrol(ofig,'style','pop','string',sprintf('%d',[0:8])','Position',bp,'Tag','SetCluster',...
        'Callback',@Update);
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(ofig,'style','pushbutton','string','plot','Position',bp,...
        'Callback',@PlotExpts);
    bp(1) = bp(1) + bp(3) + 10;

    bp(1) = 10; bp(2) = 30;
    bp(3) = cw * 8; bp(4) = ch;
    uicontrol(ofig,'style','checkbox','string','ShowN','Position',bp,'Tag','ShowN',...
        'Callback',@Update,'value',DATA.plot.shown);
    bp(2) = bp(2)+ch;
    uicontrol(ofig,'style','checkbox','string','Psych','Position',bp,'Tag','PlotPsych',...
        'Callback',@Update,'value',DATA.plot.psych);
    set(gca,'position',[0.2 0.11 0.8 0.815]);
    DATA.toplevel = ofig;

    hm = uimenu(DATA.toplevel,'Label','Extra','Tag','PlotMenu');
    uimenu(hm,'Label','&Save','Callback',@SaveAllExpt);
    uimenu(hm,'Label','&SPksXTraisl','Callback',@PlotAllSpikes);
    DATA.timerobj = timer('TimerFcn',@timerfcn, 'Period', 2.0,...
         'Tag','OnlineTimer', 'ExecutionMode','FixedSpacing','UserData',DATA.toplevel);

    DATA.lastread = 0;

    set(ofig,'UserData',DATA);
else
    DATA = get(ofig,'UserData');
end

if dattime & ~isnew
    ton = strcmp(get(DATA.timerobj,'Running'),'on');
 %don't really want timer running.  let spike2 do replots after each trial
    if ~ton
%        start(DATA.timerobj);
    end
    if ton
        stop(DATA.timerobj);
    end
    DATA.dattime = dattime;
    DATA.Expt = Expt;
    DATA.ncalls = DATA.ncalls+1;
    set(ofig,'UserData',DATA);
    PlotAllExpt(DATA.Expt);
    return;
end
if isempty(Expt) %just do plot
    PlotAllExpt(DATA.Expt);
    return;
end
DATA.Expt = Expt;
DATA.ncalls = DATA.ncalls+1;

if ~isfield(DATA,'probe')
    fprintf('Updating DATA\n');
    DATA = Update(DATA);
end

es = [Expt.Expts.Start];
if isfield(Expt.Expts,'End')
ee = [Expt.Expts.End];
else
    ee = [];
end
%would be fastest to find trials first, but need to make sure
%Expt.Trials.Start has no empties
k = 0;
for e = 1:length(expid)
for j = 1:length(Expt.Trials)
    if ~isempty(Expt.Trials(j).Start)
        starts(j) = Expt.Trials(j).Start;
        if starts(j) > es(expid(e)) && (expid(e) > length(ee) || starts(j) < ee(expid(e)));
            k = k+1;
            for p = 1:size(Expt.Trials(j).Spikes,1)
                for cl = 1:size(Expt.Trials(j).Spikes,2)
                    counts(k,p,cl) = length(Expt.Trials(j).Spikes{p,cl});
                end
            end
            Expt.Trials(j).OptionCode = '+cf';
            Expt.Trials(j).Trial = j;
            Expt.Trials(j).sz = abs(Expt.Trials(j).wi + i * Expt.Trials(j).hi);
            E.Trials(k) = Expt.Trials(j);
            sz = size(Expt.Trials(j).Spikes);
            if sz(1) >= DATA.probe & sz(2) >= DATA.cluster
            E.Trials(k).Spikes = Expt.Trials(j).Spikes{DATA.probe,DATA.cluster}' .* 10000;
            else
            E.Trials(k).Spikes = [];
            end
            np = max([sz(1) DATA.nprobes]);
            E.Trials(k).Start  = starts(j) .*10000;
            E.Trials(k).End  = Expt.Trials(j).End .*10000;
        end
    end
end
end

E.Header.Name = 'Online';
E.Header.expname = 'Online';
E.Header.rc = 0;
if isfield(E,'Trials')
E.Stimvals = Expt.Expts(expid(end)).Stimvals;
if isfield(E.Trials,'st')
    E.Stimvals.st = median([E.Trials.st]);
end
args = PlotArgs(DATA);
if sum(ismember([E.Trials.RespDir],[-1 1])) > 2
    PlotExpt(E,args{:});
    GetFigure('Psych');
    hold off;
    ExptPsych(E,'shown');
else
    PlotExpt(E,args{:});
end
%plot([E.Trials.Start],counts(:,DATA.probe,DATA.cluster),'o');
if np > DATA.nprobes
    it = findobj(DATA.toplevel,'Tag','SetProbe');
    set(it,'string',num2str([1:DATA.nprobes]'));
end
end

set(ofig,'UserData',DATA);


function timerfcn(a, varargin)
% DATA = get(findobj('Tag',get(tim,'Tag')),'UserData');
f = get(a,'UserData'); %this stores the figure
if isfigure(f)
    DATA = get(f,'UserData');
    if DATA.dattime > DATA.lastread
        GetFigure(f);
        PlotAllExpt(DATA.Expt);
        DATA.lastread = DATA.dattime;
        set(DATA.toplevel,'UserData',DATA);
    end
    stop(a);
 %fprintf('Timer Called\n');
else
    stop(a);
end


function PlotExpts(a,b)
DATA = GetDataFromFig(a);
PlotAllExpt(DATA.Expt);

function DATA = GetDataFromFig(a)

    DATA = get(a,'UserData');

    b = a;
    while isempty(DATA) & ~isempty(b)
        b = get(b,'parent');
        DATA = get(b,'UserData');
    end
    if isfield(DATA,'parentfigtag')
        DATA = get(findobj('Tag',DATA.parentfigtag),'UserData');
    end

function SaveAllExpt(a,b)
  DATA = GetDataFromFig(a);
  AllExpt = DATA.Expt;
  AllExpt.savetime = now;
  AllExpt.ncalls = DATA.ncalls;
  uisave('AllExpt','AllExpt');

function args = PlotArgs(DATA)
  args = {};
   if DATA.plot.shown
       args = {args{:} 'shown'};
   end
   if DATA.plot.psych
       args = {args{:} 'psych'};
   end
       
function DATA = Update(a,b)

if isstruct(a)
    DATA = a;
else
DATA = GetDataFromFig(a);
end
DATA.plot.shown = GetCheck('ShowN',DATA.toplevel);
DATA.probe = get(findobj(DATA.toplevel,'Tag','SetProbe'),'value');
DATA.cluster = get(findobj(DATA.toplevel,'Tag','SetCluster'),'value');
set(DATA.toplevel,'UserData',DATA);
PlotAllExpt(DATA.Expt);

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

function PlotAllSpikes(a,b)

DATA = GetDataFromFig(a);
for t = 1:length(DATA.Expt.Trials)
    for p = 1:length(DATA.Expt.Trials(t).Spikes)
        counts(t,p) = length(DATA.Expt.Trials(t).Spikes{p});
    end
end
hold off; 
imagesc(counts);
colorbar;
