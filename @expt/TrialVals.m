function X = TrialVals(E, varargin)
%X = TrialVals(Expt,varargin) Return values from trials

f = {};
j = 1;
while j <= length(varargin)
    if iscellstr(varargin{j})
        f = varargin{j};
    elseif ischar(varargin{j});
        f = {f{:} varargin{j}};
    end
    j = j+1;
end

for j = 1:length(f);
    for t = 1:length(E.Trials)
        T = E.Trials(t);
        if strcmp(f{j},'start') %nb. lower case forces it to just start(1)
            X.start(t) = T.Start(1);
        elseif strcmp(f{j},'duration')
            X.duration(t) = T.End(end)-T.Start(1);
        elseif isfield(T,f{j})
            X.(f{j})(t) = T.(f{j});
        end
    end
end