function ns = tab2im(All, varargin)

Expt = [];

j = 1;
while j < nargin
    if isfield(varargin{j},'Stimvals')
        Expt = varargin{j};
    end
    j = j+1;
end

if isinf(All(3,1))
    ns.seeds = All(1:4:end,1);
    ns.ids = All(2:4:end,1);
    ns.lsfs = All(1:4:end,2:end);
    ns.phases = All(2:4:end,2:end) * pi/180;
    ns.rsfs = All(3:4:end,2:end);
    ns.sfs = ns.lsfs & ns.rsfs;
    ns.dp = All(4:4:end,2:end) * pi/180;
    ns.dpvals = unique(ns.dp(:));
elseif isnan(All(3,1))
    ns.seeds = All(1:3:end,1);
    ns.ids = All(2:3:end,1);
    ns.sfs = All(1:3:end,2:end);
    ns.phases = All(2:3:end,2:end) * pi/180;
    ns.dp = All(3:3:end,2:end) * pi/180;
    ns.dpvals = unique(ns.dp(:));
else
    ns.seeds = All(1:2:end,1);
    ns.ids = All(2:2:end,1);
    ns.sfs = All(1:2:end,2:end);
    ns.phases = All(2:2:end,2:end) * pi/180;
end

if ~isempty(Expt)
    ns.name = splitpath(Expt.Header.Name);
    if Expt.Stimvals.st == Expt.Header.nsines
        sf = GetEval(Expt,'sf');
        df = Expt.Stimvals.f2 - sf;
        nc = Expt.Stimvals.nc;
        sfvals(1) = sf - df * nc/2;
        sfvals(nc) = sfvals(1) + df * (nc-1);
        if ismember(Expt.Stimvals.sM,[7 9 10])
            sfvals = exp(linspace(log(sfvals(1)),log(sfvals(nc)),nc));
        else
            sfvals = linspace(sfvals(1),sfvals(nc),nc);
        end
        ns.sfvals = sfvals;
    end
end