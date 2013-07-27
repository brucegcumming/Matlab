function Expt = CondenseRC(Expt, code)
%CondenseRC(Expt, code)
%takes an Expt file with RC data, i.e. many code values per trial, and
%condeses it down to 1 stimval per trial, so that PlotExpt can use it.

eb = Expt.Stimvals.e2;
ea = Expt.Stimvals.et;
Expt = FillTrials(Expt,'Dw');
if (strmatch(Expt.Stimvals.e2,'Dc') & strmatch(Expt.Stimvals.et,{'dx', 'Pd'})) ...
        | (strmatch(Expt.Stimvals.et,'Dc') & strmatch(Expt.Stimvals.e2,{'dx', 'Pd'})) ...
        | strcmp(Expt.Stimvals.e2,'e0')
    if strmatch(Expt.Stimvals.e2,{'dx', 'Pd'})
        sv = eb;
    else
        sv = ea;
    end
   for j = 1:length(Expt.Trials)
        if ~isempty(Expt.Trials(j).Start & ~isempty(Expt.Trials(j).(sv)))
            Expt.Trials(j).Start = Expt.Trials(j).Start(1);
            Expt.Trials(j).End = Expt.Trials(j).End(end);
            if ~strcmp(eb,'e0') && isfield(Expt.Trials,eb);
                Expt.Trials(j).(eb) = Expt.Trials(j).(eb)(end);
            end
            x = Expt.Trials(j).(ea)(end);
            if x > -1000
                Expt.Trials(j).(ea) = x;
            elseif x == -1006 %high x
                Expt.Trials(j).(ea) = 0;
            elseif x == -1009
                Expt.Trials(j).st = 0;
                Expt.Trials(j).(ea) = 0;
            else
                Expt.Trials(j).(ea) = 0;
                if x == -1006
                end
            end
                
            if isfield(Expt.Trials,'st') & ~isempty(Expt.Trials(j).st)
                Expt.Trials(j).st = Expt.Trials(j).st(end);
            end
            if isfield(Expt.Trials,'ce') & length(Expt.Trials(j).ce) > 1
                Expt.Trials(j).ce = Expt.Trials(j).ce(end);
            end
        else
            Expt.Trials(j).Start = NaN;
        end
    end
    rmtrials = find(isnan([Expt.Trials.Start]));
    
elseif strmatch(Expt.Stimvals.e2,'or') | strmatch(Expt.Stimvals.et,'or') 
          
for j = 1:length(Expt.Trials)
    if ~isempty(Expt.Trials(j).Start) & ~isempty(Expt.Trials(j).or)
        Expt.Trials(j).Start = Expt.Trials(j).Start(1);
        Expt.Trials(j).End = Expt.Trials(j).End(end);
        if isfield(Expt.Trials,'ori')
            Expt.Trials(j).or = Expt.Trials(j).ori;
        elseif(Expt.Trials(j).Dw < 18)
            ors = unique(Expt.Trials(j).or(1:end-1));
            Expt.Trials(j).or = Expt.Trials(j).or(end);
            exptor(j) = Expt.Trials(j).or;
        end
    else
        Expt.Trials(j).Start = NaN;
    end
    if isfield(Expt.Trials,'me')
        Expt.Trials(j).me = median(Expt.Trials(j).me);
    end
    if isfield(Expt.Trials,'st')
        Expt.Trials(j).st = median(Expt.Trials(j).st);
    end
end
if Expt.Trials(1).Dw >=18
    Expt.Trials(1).or = mode(exptor);
end
rmtrials = find(isnan([Expt.Trials.Start]));
for dw = unique([Expt.Trials.Dw])
    id = find([Expt.Trials.Dw] == dw);
    if length(id) < length(Expt.Trials)/20
        rmtrials = [rmtrials id];
    end
end
for j = 2:length(Expt.Trials)
    if(Expt.Trials(j).Dw >= 18 & ~isempty(Expt.Trials(j).or))
        Expt.Trials(j).or = Expt.Trials(j-1).or;
    end
end
else
    sv = Expt.Stimvals.et;
   for j = 1:length(Expt.Trials)
        if ~isempty(Expt.Trials(j).Start & ~isempty(Expt.Trials(j).(sv)))
            Expt.Trials(j).Start = Expt.Trials(j).Start(1);
            Expt.Trials(j).End = Expt.Trials(j).End(end);
            Expt.Trials(j).(eb) = Expt.Trials(j).(eb)(end);
            Expt.Trials(j).(ea) = Expt.Trials(j).(ea)(end);
            if isfield(Expt.Trials,'st')
                Expt.Trials(j).st = Expt.Trials(j).st(end);
            end
            if isfield(Expt.Trials,'ce')
                Expt.Trials(j).ce = Expt.Trials(j).ce(end);
            end
            if isfield(Expt.Trials,'me')
                Expt.Trials(j).me = Expt.Trials(j).me(end);
            end
        else
            Expt.Trials(j).Start = NaN;
        end
    end
    rmtrials = find(isnan([Expt.Trials.Start]));
    
    
end


if ~isempty(rmtrials)
    id = setdiff([1:length(Expt.Trials)],rmtrials);
    Expt.Trials = Expt.Trials(id);
end




