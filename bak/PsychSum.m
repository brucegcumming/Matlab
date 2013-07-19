function ps = PsychSum(name, varargin)
%ps = PsychSum(name, varargin)
%read all online psych files from directory name
%and calculate psfs. By default only does NoiseOri files
%
%PsychSum(name, 'prefix', prefix)
%uses files starting with prefix in the name

mintrials = 100;
et = 'or';
sincedate = 0;
rebuild = 0;
saveps = 0;
outfile = [];
prefix = 'NoiseOri'; 
plotargs = {};
j = 1;
args = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'absmax',6)
        args = {args{:} 'xmax' abs(varargin{j}) 'xmin' -abs(varargin{j})};
    elseif strncmpi(varargin{j},'hold',4)
        plotargs = {plorargs{:} 'hold'};
    elseif strncmpi(varargin{j},'prefix',5)
        j = j+1;
        prefix = varargin{j};
        if strncmpi(prefix,'Plaid',4)
            et = 'pR';
        end
    elseif strncmpi(varargin{j},'rebuild',5)
        rebuild = 1;
    elseif strncmpi(varargin{j},'save',4)
        saveps = 1;
    elseif strncmpi(varargin{j},'since',5)
        j = j+1;
        sincedate = datenum(varargin{j});
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end

np = 1;
if iscellstr(name)
    names = name;
elseif isstruct(name)
    PlotThs(name, varargin{:});
    return;
elseif isdir(name)
        [names, sizes, dates] = TreeFind(name,'name',prefix);
end


if rebuild == 0
for j = 1:length(names)
    psum = [];
    if strfind(names{j},'.mat') & strfind(names{j},prefix)
        load(names{j});
        sincedate = max([psum.StartTime])+0.5;
        if sincedate > now
            dates = [psum.StartTime];
            sincedate = max(dates(dates < now));
        end
        outfile = names{j};
    end
end
end

if sincedate
    id = find(dates > sincedate);
    names = names(id);
end
good = [];
for j = 1:length(names)
      if isempty(regexp(names{j},'.*.log')) && isempty(regexp(names{j},'.*.mat'))
         good = [good j];
      end
end
names = names(good);
ps = [];
for j = 1:length(names)
        exs = PsychMon(names{j},'getexpts');
        types = [];
        btypeid = [];
        btypes = {};
        for e = 1:length(exs)
            Expt=exs{e}; 
            if strcmp(Expt.Stimvals.et,et)
            if length(Expt.Trials) > mintrials
                types(e) = max([Expt.Trials.(et)]);
            else
                types(e) = 0;
            end
            else
                types(e) = NaN;
            end
            btypes{e} = Expt.Stimvals.e2;
        end
        st = unique(types);
        st = st(st > 0);
        eb = unique(btypes);
        for e = 1:length(exs)
            btypeid(e) = strmatch(btypes{e},eb);
        end
        n = 0;
        allst = [];
        alleb = [];
        for e = 1:length(st)
            for b = 1:length(eb)
                id = find(types == st(e) & btypeid == b);
                if length(id)
                    n = n+1;
                    allst(n) = st(e);
                    alleb(n) = strmatch(eb{b},eb);
                end
            end
        end
        st = allst;
        eb = alleb;
        fprintf('%s %d/%d expts',names{j},length(st),e);
        fprintf('%.1f ,',st);
        fprintf('\n');
        for e = 1:length(st)
            id = find(types == st(e) & btypeid == eb(e));
            GetFigure('Psych');
            hold off;
            [pp, p] = ExptPsych(exs(id),'sum','shown', args{:});
            drawnow;
            if length(pp) > 1
            if isfield(exs{id(1)}.Stimvals,'Start')
                p.StartTime = exs{id(1)}.Stimvals.Start;
            else
                p.StartTime = NaN;
            end
            f = {'fit' 'maxperf' 'labels' 'rws' 'rwn' 'StartTime'};
            for fi = 1:length(f)
            ps(np).(f{fi}) = p.(f{fi});
            end
            ps(np).e1val = st(e);
            ps(np).e2type = btypes{id(1)};
            ps(np).filename = names{j};
            thr(np) = ps(np).fit.fit(1);
            np = np+1;
            end
        end
end
if rebuild == 0
    ps = [psum ps];
end
PlotThs(ps, plotargs{:});
if saveps && length(outfile)
    fprintf('Saving %s\n',outfile);
    psum = ps;
    save(outfile,'psum');
end

function HitPoint(a,b, id)

DATA = GetDataFromFig(a);
fprintf('%d %s\n',id,DATA.ps(id).filename);
GetFigure('PSF');
if DATA.holdon
    hold on;
else
hold off;
end
fitpsf(DATA.ps(id).fit.data,'showfit',DATA.ps(id).fit,'shown');
title(sprintf('%s %.2f %s',DATA.ps(id).labels{:},abs(DATA.ps(id).fit.fit(2)),datestr(DATA.ps(id).StartTime)));




function PlotThs(ps, varargin)

plottype = 'Default';

id = 1:length(ps);
c = struct2cell(ps);
DATA.holdon = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'hold',4)
        DATA.holdon = 1;
    elseif strncmpi(varargin{j},'since',5)
        j = j+1;
        sincedate = datenum(varargin{j});
        id = find([ps.StartTime] > sincedate);
    elseif strncmpi(varargin{j},'Type2',5)
        j = j+1;
        tid = strmatch(varargin{j},c(8,1,id));
        id = id(tid);
    elseif strncmpi(varargin{j},'TimeCourse',4)
        plottype = 'TimeCourse';
    end
    j = j+1;
end

for e = 1:length(ps)
            thr(e) = abs(ps(e).fit.fit(2));
            ors(e) = max([ps(e).fit.data.y]);
            starts(e) = ps(e).StartTime;
end
GetFigure('Psych Summary');
if strmatch(plottype,'TimeCourse')
    subplot(1,1,1);
    hold off;
    plot(starts(id),ors(id),'o');
    hold on; 
    for j = 1:length(id)
        plot(starts(id(j)),thr(id(j)).*100,'ro','buttondownfcn',{@HitPoint, id(j)});
    end
    datetick('x','dd/mm');
    set(gca,'yscale','log');
else
subplot(2,1,1);
plot(ors(id),thr(id),'o');
subplot(2,1,2);
hold off;
plot(starts(id),ors(id),'o');
hold on; 
for j = 1:length(id)
    plot(starts(id(j)),thr(id(j)).*100,'ro','buttondownfcn',{@HitPoint, id(j)});
end
datetick('x','dd/mm');
set(gca,'yscale','log');
end
DATA.ps = ps;
DATA.thr = thr;
DATA.ors = ors;
DATA.starts = starts;
set(gcf,'UserData',DATA);