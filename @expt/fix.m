function Expt = fix(Expt, type, varargin)
% expt.fix inconsistencies in Expt structs
%Expt = expt.fix(Expt,'starts') ensures start times are monotonically
%increasing
% 

if nargin == 1
    type = 'all';
end

if iscell(Expt)
    for j = 1:length(Expt)
        Expt{j} = expt.fix(Expt{j}, type, varargin{:});
    end
    return;
end

if sum(strcmp(type,{'all' 'starts'}))
    toff = 0;
    last = 0;
    for t = 1:length(Expt.Trials)
        if Expt.Trials(t).Start(1)+toff < last
            toff = 1000 + last - Expt.Trials(t).Start(1); %last has current toff in it
        end
        Expt.Trials(t).Start = Expt.Trials(t).Start + toff;
        last =  Expt.Trials(t).Start(end);
        if toff > 0
            Expt.Trials(t).timeoffset = toff;
        end
    end
end