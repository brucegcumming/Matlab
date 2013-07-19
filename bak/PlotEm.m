function PlotEm(Expt, trials, varargin)
%
%PlotEm(Expt, trials, ......)
%Plots eye movement traces from .em mat files
%
% N.B. This is different from the stucture of Expt files to which eye
% movement data have been added - those need PlotExptEM(....
%
%  Expt is the Expt structure containing Trials and Header etc.
%  trials is a list of trial ids
%
% PlotEm(Expt, 1)  Plots the data from the first trial
% PlotEm(Expt, 1:3)  Plots the data from the first three trials,
% superimposed
% the default plot is 4 traces wrt time - lh,rh,lv,rv
%
% PlotEm(Expt, 1, 'xy') plots v,h positions on an X-Y plot
%  PlotEm(Expt, 1, 'cv') conjugate H, H vergence, and conjugate V
%

plottype = 0;
plotavg = 0;
plotmeans = 0;
plotsac = 0;
figlabel{1} = 'EmPlotYT';
showtime = 0;
nohold = 1;
eyes = 2;
yrange = [];
Text = [];
calib = [1 1 1 1 1 1; 0 0 0 0 0 0];

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'text')
            Text = varargin{j};
        end
    elseif strncmpi(varargin{j},'rawavg',3)
        plotavg = 1;
    elseif strncmpi(varargin{j},'avg',3)
        plotavg = 2;
    elseif strncmpi(varargin{j},'calib',5)
        j = j+1;
        calib = varargin{j};
    elseif strncmpi(varargin{j},'cv',2)
        plottype = 1;
    elseif strncmpi(varargin{j},'hold',2)
        nohold = 0;
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        ids = varargin{j};
        trials = find(ismember([Expt.Trials.id],ids));
    elseif strncmpi(varargin{j},'means',4)
        plotmeans = 1;
    elseif strncmpi(varargin{j},'three',4)
        eyes = 3;
    elseif strncmpi(varargin{j},'speed',2)
        plottype = 3;
    elseif strncmpi(varargin{j},'yrange',4)
        j = j+1;
        yrange = varargin{j};
    elseif strncmpi(varargin{j},'saccades',3)
        plotsac = 1;
    elseif strncmpi(varargin{j},'xy',2)
        plottype = 2;
    end
    j = j + 1;
end

[efig, isnew] = GetFigure(figlabel{1});
if isnew
    DATA.trial = trials(end);
    DATA.Expt = Expt;
    DATA.Text = Text;
    DATA.yrange = yrange;
    DATA.calib = calib;
    DATA.neyes = eyes;
        bp = [10 10 40 20];
    uicontrol(efig,'style','pushbutton','string','>>','Position',bp,'Tag','NextTrial',...
        'Callback', {@PlayNextTrial, 1});
    bp(1) = bp(1)+bp(3);
    uicontrol(efig,'style','pushbutton','string','<<','Position',bp,'Tag','PrevTrial',...
        'Callback', {@PlayNextTrial, -1});
    set(gcf,'UserData',DATA);
else  %may be setting something
    DATA = get(efig,'UserData');
    if ~isempty(yrange)
        DATA.yrange = yrange;
    end
    set(efig,'UserData',DATA);
end

if isempty(trials)
    return;
end
if ~isempty(Text) & plotsac == 1
    PlotTrialSaccades(DATA, trials);
    return;
end

if ~isfield(Expt.Header,'CRsamplerate')
    Expt.Header.CRsamplerate = 1./(Expt.Header.CRrates(1) * 10000);
end
if ~isfield(Expt.Trials,'Trial')
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).Trial = j;
    end
end
scales = DATA.calib(1,:);

if plotavg
    buflen = NaN;
    for j = 1:length(trials);
        buflen = min([buflen length(Expt.Trials(j).Eyevals.lh)]);
    end
    lvs = [];
    rvs = [];
    lhs = [];
    rhs = [];
    for j = 1:length(trials);
        trial = trials(j);
        lvs(j,:) = Expt.Trials(trial).Eyevals.lv(1:buflen)';
        rvs(j,:) = Expt.Trials(trial).Eyevals.rv(1:buflen)';
        rhs(j,:) = Expt.Trials(trial).Eyevals.rh(1:buflen)';
        lhs(j,:) = Expt.Trials(trial).Eyevals.lh(1:buflen)';
        if plotavg == 2
            lvs(j,:) = lvs(j,:) - mean(lvs(j,:));
            rvs(j,:) = rvs(j,:) - mean(rvs(j,:));
            rhs(j,:) = rhs(j,:) - mean(rhs(j,:));
            lhs(j,:) = lhs(j,:) - mean(lhs(j,:));
        end
    end
    lvs(find(lvs < -99)) = NaN;
    rvs(find(rvs < -99)) = NaN;
    ev.lh = mean(lhs)';
    ev.rh = mean(rhs)';
    ev.lv = mean(lvs)';
    ev.rv = mean(rvs)';
    ntrials = 1;
elseif plotmeans
    times = [1:3000]./(Expt.Header.CRsamplerate * 10000); %% Time in sec, not tics
    ts = find(times < 2.0);
    for j = 1:length(Expt.Trials)
        if length(Expt.Trials(j).Eyevals.rh) >= ts(end)
        h(j) = mean(Expt.Trials(j).Eyevals.rh(ts));
        v(j) = mean(Expt.Trials(j).Eyevals.rv(ts));
        end
        good(j) = Expt.Trials(j).good;
        hs(j) = mean(Expt.Trials(j).softoff(1));
        vs(j) = mean(Expt.Trials(j).softoff(3));
    end
    hold off;
    id = find(good ==1);
    plot(h(id),'o');
    hold on;
    plot(v(id),'ro');
    plot(hs(id),'+');
    plot(vs(id),'r+');
else
    ntrials = length(trials);
end


reclen = 0;
if nohold
    hold off;
else
    hold on;
end
for j = 1:ntrials;
    if ~plotavg
        if isfield(Expt.Trials,'EyeData')
            ev.lh = Expt.Trials(trials(j)).EyeData(:,1);
            ev.rh = Expt.Trials(trials(j)).EyeData(:,2);
            ev.lv = Expt.Trials(trials(j)).EyeData(:,3);
            ev.rv = Expt.Trials(trials(j)).EyeData(:,4);

        else
            ev = Expt.Trials(trials(j)).Eyevals;
        end
    end
% sometimes # samples differs by 1
    reclen = min([length(ev.lh) length(ev.lv) length(ev.rh) length(ev.rv)]);
    if eyes == 3
    rid = 1:(reclen-1);
    else
        rid = 1:reclen;
    end
    times = rid./(Expt.Header.CRsamplerate * 10000); %% Time in sec, not tics
    if plottype == 1
        plot(times,ev.lh-ev.rh,'r');
        hold on;
        plot(times,mean([ev.lh ev.rh],2),'g');
        plot(times,mean([ev.lv ev.rv],2),'b');
        if plotavg
            plot(times,ev.lh-ev.rh+std(lhs - rhs)','r:');
        end
    elseif plottype == 2
        if showtime
            for k = 1:length(ev.lh);
                plot(ev.lh(1:k),ev.lv(1:k),'r');
                hold on;
                plot(ev.rh(1:k),ev.rv(1:k),'r');
                if showtime & mod(k,showtime) == showtime-1
                    drawnow;
                end
            end
        else
            plot(ev.lh,ev.lv,'r');
            hold on;
            plot(ev.rh,ev.rv,'g');
        end
    elseif plottype == 3
        if eyes == 3
            hold off;
            plot(eyespeed(ev,2));
            hold on;
            plot(eyespeed(ev,3),'r');
        else
            plot(eyespeed(ev,0));
        end
    else
        plot(times,ev.lh(rid).*scales(1),'r');
        hold on;
        plot(times,ev.rh(rid).*scales(2),'g');
        plot(times,ev.lv(rid).*scales(3),'m');
        plot(times,ev.rv(rid).*scales(4),'c');
        if eyes == 3
            plot(times,ev.xv(rid).*scales(6),'c');
            plot(times,ev.xh(rid).*scales(5),'g');
        end
    end
end
DATA.trial = trials(1);
if plottype == 0
    xlabel('Time (sec)');
    legend('LH','RH','LV','RV');
elseif plottype == 1
    xlabel('Time (sec)');
    legend('H','V','H verg');
end
if ~isempty(yrange)
    set(gca,'ylim',yrange);
end
if length(trials) < 10
    title(['Trials  ' sprintf('%d,',[Expt.Trials(trials).Trial])]);
else
    title(sprintf('%d Trials %d-%d ',length(trials),min([Expt.Trials(trials).Trial]),max([Expt.Trials(trials).Trial])));
end
set(gcf,'UserData',DATA);


function PlayNextTrial(a,b,step)

DATA = get(gcf,'UserData');
if DATA.trial+step <= length(DATA.Expt.Trials) && DATA.trial+step > 0
DATA.trial = DATA.trial+step;
end
hold off;
PlotTrialSaccades(DATA, DATA.trial);
if isfield(DATA.Expt.Trials,'fx') && isfield(DATA.Expt.Trials,'fy')
    xs = sprintf('(%.1f,%.1f)',DATA.Expt.Trials(DATA.trial).fx,DATA.Expt.Trials(DATA.trial).fy);
else
xs='';
end
title(sprintf('Trial %d (id%d) at %.3f%s',DATA.trial,DATA.Expt.Trials(DATA.trial).id,DATA.Expt.Trials(DATA.trial).Start(1)/10000,xs));
set(gcf,'UserData',DATA);

function speed = eyespeed(ev, eye)

if eye == 0
ch = diff(smooth((ev.lh + ev.rh)/2,2));
cv = diff(smooth((ev.lv + ev.rv)/2,2));
elseif eye == 2
ch = diff(smooth((ev.rh),2));
cv = diff(smooth((ev.rv),2));
elseif eye == 3
ch = diff(smooth((ev.xh),2));
cv = diff(smooth((ev.xv),2));
end
    speed = sqrt(ch.^2 + cv.^2);
   
    
function PlotTrialSaccades(DATA, trials, varargin)

Expt = DATA.Expt;
if isfield(DATA,'Text') && ~isempty(DATA.Text)
Text = DATA.Text;
sid = strmatch('Sa',Text.text);
else
    sid = [];
    id = [];
end
scales = DATA.calib(1,:);
offsets = DATA.calib(2,:);
trials = trials(find(trials < length(Expt.Trials)));
for j = 1:length(trials)
    t = trials(j);
    if isfield(Expt.Trials,'ftime')
        pre = Expt.Trials(t).Start(1)/10000 - Expt.Trials(t).ftime/10000;
        start = Expt.Trials(t).ftime./10000;
    else
        pre = 0.1;
        start = Expt.Trials(t).Start(1)./10000 - pre;
    end
    if isfield(Expt.Trials,'EyeData')
        Eyevals.lh = Expt.Trials(trials(j)).EyeData(:,1) .* scales(1);
        Eyevals.rh = Expt.Trials(trials(j)).EyeData(:,2) .* scales(2);
        Eyevals.lv = Expt.Trials(trials(j)).EyeData(:,3) .* scales(3);
        Eyevals.rv = Expt.Trials(trials(j)).EyeData(:,4) .* scales(4);
        npts = size(Expt.Trials(trials(j)).EyeData,1);
    else
        npts = min([length(Expt.Trials(t).Eyevals.rh) length(Expt.Trials(t).Eyevals.lh) length(Expt.Trials(t).Eyevals.rv) length(Expt.Trials(t).Eyevals.lv)]);
        Eyevals = Expt.Trials(t).Eyevals;
    end
    ts = ([1:npts] * Expt.Header.CRrates(1))-pre;
    smps = 1:npts;
    endtime = start+max(ts);
    plot(ts,Eyevals.rh(smps).*scales(2),'g');
    hold on;
    plot(ts,Eyevals.lh(smps).*scales(1),'r');
    plot(ts,Eyevals.lv(smps).*scales(3),'m');
    plot(ts,Eyevals.rv(smps).*scales(4),'c');
    if DATA.neyes == 3
    plot(ts,Eyevals.xh(smps).*scales(5),'g');
    plot(ts,Eyevals.xv(smps).*scales(6),'m');
    end
    
    if length(sid)
    id = find(Text.times(sid)> start & Text.times(sid)< endtime);
    ts = ([1:length(Eyevals.rv)] * Expt.Header.CRrates(1))-pre;
    plot(ts,Eyevals.rv,'r');
    end
    if length(DATA.yrange) == 2
        set(gca,'ylim', DATA.yrange);
        yl = DATA.yrange;
    else
        yl = get(gca,'ylim');
    end
    if isfield(Expt.Trials,'Saccades')
        tr = trials(j);
        for k = 1:length(Expt.Trials(tr).Saccades)
            [a, ti] = min(abs(ts-(Expt.Trials(tr).Saccades(k).peakt./10000-pre)));
            y(1) = 0;
            y(2) = Expt.Trials(tr).Saccades(k).size;
            plot([ts(ti) ts(ti)],y,'r')
        end
    end
    for k = 1:length(id)
        ss =sscanf(Text.text(sid(id(k)),:),'Sa%f %f');
        if length(ss) == 2
        ti = round((Text.times(sid(id(k)))-start)./Expt.Header.CRrates(1));
        y(1) = Eyevals.rh(ti);
        y(2) = Eyevals.rh(ti)-ss(1);
        plot([ts(ti) ts(ti)],y,'b')
        y(1) = Eyevals.rv(ti);
        y(2) = Eyevals.rv(ti)-ss(2);
        plot([ts(ti) ts(ti)],y,'r')
        h = text(ts(ti),yl(2),sprintf('%.2f\n%.2f',ss(1),ss(2)));
        set(h,'VerticalAlignment','Top');
        fprintf('Saccade %.2f,%.2f at %.3f\n',ss(1),ss(2),ts(ti)+start);
        end
    end
    if isfield(DATA,'Text') && ~isempty(DATA.Text)

    id = find(Text.times > start & Text.times < endtime & Text.codes(:,4) == 2);
    for k = 1:length(id)
        ti = Text.times(id(k))-start;
        if Text.codes(id(k),1) == 11  %Badfix
            plot([ti ti],yl,'c');
        elseif Text.codes(id(k),1) == 4  %Badfix
            plot([ti ti],yl,'m');
        end
    end
    end
end
   
function PlotSaccades(Expt, Text, trials, varargin)
        
starts = [Expt.Trials.ftime]./10000;
ends = [Expt.Trials.End]./10000;
sid = strmatch('Sa',Text.text);
for j = 1:length(trials)
    s = sid(trials(j));
    t = Text.times(s);
    id = find(starts <t & ends > t);
    if length(id)
        ss =sscanf(Text.text(s,:),'Sa%f %f');
        k = id(1);
        ts = round((t-starts(k))/Expt.Header.CRrates(1));
        if ts > 100 & ts < length(Expt.Trials(k).Eyevals.rh)-100
        plot(Expt.Trials(id(1)).Eyevals.rh(ts-100:ts+100));
        hold on;
        y(1) = Expt.Trials(id(1)).Eyevals.rh(ts);
        y(2) = Expt.Trials(id(1)).Eyevals.rh(ts)-ss(1);
        plot([100 100],y,'b')
        plot(Expt.Trials(id(1)).Eyevals.rv(ts-100:ts+100),'r');
        y(1) = Expt.Trials(id(1)).Eyevals.rv(ts);
        y(2) = Expt.Trials(id(1)).Eyevals.rv(ts)-ss(2);
        plot([100 100],y,'r');
        end
    end
end