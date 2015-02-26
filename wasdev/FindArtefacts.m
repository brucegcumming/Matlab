function result = FindArtefacts(DATA, varargin)

stimes= [];
first=1;
second = 2;
mint = 0.00025;
findcount = [];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'count',5)
        j = j+1;
        findcount = varargin{j};
    elseif strncmpi(varargin{j},'mint',4)
        j = j+1;
        mint = varargin{j};
    end
    j = j+1;
end

if isstruct(DATA)
    result = PlotResults(DATA, varargin{:});
    return;
end

if ischar(DATA) 
    spkdir = [DATA '/Spikes'];
    [a, num] = fileparts(DATA);
    DATA = [];
    [b, monk] = fileparts(a);
    if exist(spkdir,'dir')
        d = dir(spkdir);
        prefs = {};
        for j = 1:length(d)
            pref = regexprep(d(j).name,'t[0-9]\..*','');
            if strncmp(pref,monk,3) & isempty(strmatch(pref,prefs))
                prefs = {prefs{:} pref};
            end
        end
        t = 9;
        for j = 1:length(prefs);
            sfile{j} = sprintf('%s/%st%d.mat',spkdir,prefs{j},t);
            if exist(sfile{j},'file')
                A = load(sfile{j});
                f = fields(A);
                id = strmatch('Ch',f);
                vars{j} = f{id};
                p = sscanf(A.(f{id}).title,'Spike %d');
                DATA.AllSpikes{p} = A.(f{id});
                probes(j) = p;
                clear A;
            end
        end
    end
end

tstart = now;
for lap = 1:12
    fprintf('Lap %d\n',lap);
    first = lap*2-1;
    second = lap*2;
for j = 1:length(DATA.AllSpikes)
sid{j} = [];
end
for j = 1:length(DATA.AllSpikes{first}.times)
    id = find(abs(DATA.AllSpikes{second}.times-DATA.AllSpikes{first}.times(j)) < mint);
    if length(id)
        stimes = [stimes DATA.AllSpikes{first}.times(j)];
        alltimes(first,length(stimes)) = j;
        alltimes(second, length(stimes)) = id(1);
        allmeans(first,length(stimes)) = mean(DATA.AllSpikes{first}.values(j,:));
        allmeans(second,length(stimes)) = mean(DATA.AllSpikes{second}.values(id(1),:));
    end
end
fprintf('%d events',length(stimes));
counts = ones(size(stimes)).*2; %already have two
others = setdiff([1:length(DATA.AllSpikes)],[first second]);
for j = others
    fprintf('%dx%d times',size(alltimes,1),size(alltimes,2));
    fprintf('%s probe %d\n',datestr(now),j);
    for k = 1:length(stimes)
        id = find(abs(DATA.AllSpikes{j}.times-stimes(k)) < mint);
        if length(id)
        counts(k) = counts(k) +1;
        alltimes(j,k) = id(1);
        allmeans(j,k) = mean(DATA.AllSpikes{j}.values(id(1),:));
        end
    end
end
end

colors = mycolors;
for j = 1:length(stimes)
    if counts(j) > 15
    for p = first:size(alltimes,1)
        if alltimes(p,j) > 0
        plot(DATA.AllSpikes{p}.values(alltimes(p,j),:)+p*5,'color',colors{1+mod(p-1,4)});
        hold on;
        end
    end
    end
end

for j = 1:length(stimes)
    if counts(j) > 12
        for p = 1:size(alltimes,1)
            if alltimes(p,j) > 0
                DATA.AllSpikes{p}.times(alltimes(p,j)) = -abs(DATA.AllSpikes{p}.times(alltimes(p,j)));
            end
        end
    end
end
for j = 1:length(DATA.AllSpikes)
    eval([vars{j} '=DATA.AllSpikes{j};']);
    fprintf('Writting %s to %s\n',vars{j},sfile{j});
    save(sfile{j},vars{j});
end

    
result.counts = counts;
result.times = stimes;
result.means = allmeans;
result.ids = alltimes;
result.AllSpikes = DATA.AllSpikes;        
    

function result = PlotResults(res, varargin)

ids = [];
findcount = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'count',5)
        j = j+1;
        findcount = varargin{j};
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        ids = varargin{j};
    end
    j = j+1;
end
   

if isempty(ids)
    if ~isempty(findcount)
        ids = find(ismember(res.counts,findcount));
    else
    ids = find(res.counts > 16);
    end
end

[f, isnew] = GetFigure('Artifact Waveforms');
if isnew
    uicontrol(f, 'style', 'pushbutton', 'string', '>>', 'units', 'normalized',...
        'position', [0 0 0.05 0.05],'callback',{@NextWave, 1});
    uicontrol(f, 'style', 'pushbutton', 'string', '<<', 'units', 'normalized',...
        'position', [0.05 0 0.05 0.05],'callback',{@NextWave, -1});
    uicontrol(f, 'style', 'pushbutton', 'string', 'clr', 'units', 'normalized',...
        'position', [0.1 0 0.05 0.05],'callback',{@NextWave, 0});
end

res.selid = ids;
res.current = 1;
set(f,'UserData',res);
for j = 1:length(ids)
    hold off;
    result.energy(j) = PlotWave(res,ids(j));
    drawnow;
end

function e = PlotWave(res, t)
colors = mycolors;
    id = find(res.ids(:,t) > 0);
    for p = 1:length(id)
        plot(res.AllSpikes{id(p)}.values(res.ids(id(p),t),:),'color',colors{id(p)});
        mw(p,:) = res.AllSpikes{id(p)}.values(res.ids(id(p),t),:);
        hold on;
        title(sprintf('%d at %.2f (%s)',res.counts(t),res.times(t),sprintf('%d ',id)));
    end
    plot(mean(mw),'k','linewidth',2);
    e = sum(diff(mean(mw)).^2);

function NextWave(a,b, step)

f = get(a,'parent');
res = get(f,'UserData');
res.current = res.current+step;
if step == 0
    hold off;
end
PlotWave(res, res.selid(res.current));
set(f,'UserData',res);
