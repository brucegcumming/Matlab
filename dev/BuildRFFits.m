function fits = BuildRFFits(dirname)


if iscellstr(dirname)
    for j = 1:length(dirname)
        fits{j} =BuildRFFits(dirname{j});
    end
    return;
elseif iscell(dirname)
    fits = CombineFits(dirname);
    return;
end
d = mydir([dirname '/*PP.mat']);
fits = FitExpts(d);
d = mydir([dirname '/*OP.mat']);
a = FitExpts(d);
fits = {fits{:} a{:}};
fits = CombineFits(fits);

function fits = CombineFits(fits)
probe = [];
for j = 1:length(fits)
    if isfield(fits{j},'pvar')
        pvar(j) = fits{j}.pvar;
        probe(j) = fits{j}.probe;
        meanpos(j) = fits{j}.mean;
        expts{j} = fits{j}.expt;
        depth(j) = fits{j}.depth;
        if strfind(expts{j},'Pp')
            exptype(j) = 1;
        elseif strfind(expts{j},'Op')
            exptype(j) = 2;
        else
            exptype(j) = 0;
        end
        ros(j) = fits{j}.Ro;
    else
        pvar(j) = NaN;
    end
end
if ~isempty(probe)
    probes = unique(probe);
    nrf = 0;
    for j = 1:length(probes)
        pid = find(probe == probes(j) & exptype == 1);
        oid = find(probe == probes(j) & exptype == 2);
        if ~isempty(oid) && ~isempty(pid)
            d = sort(depth([oid pid]));
            bid = find(diff((d)) > 0.5); %500uM break in position
            bid = [bid length(d)];
            start = 0;
            for nd = 1:length(bid)
                dpid = pid(find(depth(pid) > start & depth(pid)  <= d(bid(nd))));
                doid = oid(find(depth(oid) > start & depth(oid)  <= d(bid(nd))));
                ppos = mean(meanpos(dpid));
                opos = mean(meanpos(doid));
                rf = op2xy([opos ppos],ros(dpid(1)));
                nrf = nrf+1;
                proberf(nrf,:) = [probes(j) rf];
                fits{end}.proberf = proberf;
                start = d(bid(nd));
            end
        end
    end
end

function fits = FitExpts(d)
fits = {};
nf = 0;
for j = 1:length(d)
    if isempty(strfind(d(j).name,'.lfp.'))
        res = PlotExpt(d(j).name,'noplot');
        nf = nf+1;
        if ~isempty(res)
        E = res.Data;
        fits{nf} = FitExpt(res);
        if isfield(E.Header,'probe')
            p = E.Header.probe;
        else
            p = GetProbeFromName(d(j).name);
        end
        fits{nf}.probe = p;
        fits{nf}.name = d(j).name;
                
        if isfield(E.Header,'depth')

            depth = E.Header.depth;
        elseif isfield(E.Stimvals,'ed')
            depth = E.Stimvals.ed;
        elseif isfield(E.Trials,'ed')
            depth = mean([E.Trials.ed]);
        else
            depth = NaN;
        end
        fits{nf}.depth = depth;
        fits{nf}.expt = E.Header.expname;
        fits{nf}.Ro = E.Stimvals.Ro;
        end
    end
end

