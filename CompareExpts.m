function CompareExpts(Ea,Eb, varargin)
%CompareExpts(Ea,Eb, varargin) Compares Expt struct and reports differences
%Fields that are always different are not reported 'all' is used.
%'all' is useful for comparing two Expts that should be identical
%CompareExpts(....'header') Does not compare indivudal trials

headeronly = 0;

j = 1;
IgnoreTrialFields = {'rptframes' 'Trial' 'End' 'RespDir' 'sq'};
IgnoreStimvals = {'TrueEnd' 'rwdir' };
while j <= length(varargin)
    if strncmpi(varargin{j},'all',3)
        IgnoreTrialFields = {};
        IgnoreStimvals = {};
    elseif strncmpi(varargin{j},'Header',3)
        headeronly = 1;
    elseif strncmpi(varargin{j},'matchtimes',7)
        MatchExptTimes(Ea, Eb);
        return;
    end
    j = j+1;
end

if iscell(Ea) && iscell(Eb)
    for j = 1:length(Ea)
        CompareExpts(Ea{j},Eb{j});
    end
    return;
end

f = setdiff(fields(Ea.Stimvals),IgnoreStimvals);
for j = 1:length(f);
    if isfield(Eb.Stimvals,f{j})
    if Ea.Stimvals.(f{j}) ~= Eb.Stimvals.(f{j})
        if sum(~isnan(Ea.Stimvals.(f{j}))) || sum(~isnan(Eb.Stimvals.(f{j}))) 
            fprintf('Stimval %s different %.2f %.2f \n',f{j},Ea.Stimvals.(f{j}),Eb.Stimvals.(f{j}));
        end
    end
    else
       fprintf('No %s in Expt2.Stimvals\n',f{j});
    end
end

f = setdiff(fields(Ea.Trials),fields(Eb.Trials));
for j = 1:length(f)
    fprintf('Trial.%s not in Expt2\n',f{j});
end
f = setdiff(fields(Eb.Trials),fields(Ea.Trials));
for j = 1:length(f)
    fprintf('Trial.%s not in Expt1\n',f{j});
end

if headeronly
    return;
end

f = setdiff(fields(Ea.Trials),IgnoreTrialFields);

for j = 1:length(f);
    if isfield(Eb.Trials,f{j})
        for t = 1:length(Ea.Trials)
            if iscell(Ea.Trials(t).(f{j}))
            elseif length(Ea.Trials(t).(f{j})) ~= length(Eb.Trials(t).(f{j}))
                fprintf('Trial %d length(%s) %d vs %d\n',t,f{j},length(Ea.Trials(t).(f{j})),length(Eb.Trials(t).(f{j})));
            elseif Ea.Trials(t).(f{j}) ~= Eb.Trials(t).(f{j})
                fprintf('Trial %d.%s %.4f vs %.4f\n',t,f{j},Ea.Trials(t).(f{j}),Eb.Trials(t).(f{j}));
            end
        end
    else
        fprintf('Expt2.Trials no %s\n',f{j});
    end
end


function MatchExptTimes(Ea, Eb, varargin)


for j = 1:length(Ea.Trials)
    Sa(j) = Ea.Trials(j).Start(1);
end
for j = 1:length(Eb.Trials)
    Sb(j) = Eb.Trials(j).Start(1);
end
Ia = [Ea.Trials.id];
Ib = [Eb.Trials.id];

hold off;
plot(Sa, [Ea.Trials.id],'o');
hold on;
plot(Sb, [Eb.Trials.id],'ro');
for j = 1:length(Sa)
    id = find(Sb == Sa(j) & Ib == Ia(j)) ;
    if isempty(id)
        plot(Sa(j),Ia(j),'bo','markerfacecolor','b');
    end
end
for j = 1:length(Sb)
    id = find(Sa == Sb(j) & Ia == Ib(j)) ;
    if isempty(id)
        plot(Sb(j),Ib(j),'ro','markerfacecolor','r');
    end
end
