function CompareExpts(Ea,Eb, varargin)

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'matchtimes',7)
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

f = fields(Ea.Stimvals);
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
f = fields(Ea.Trials);
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
