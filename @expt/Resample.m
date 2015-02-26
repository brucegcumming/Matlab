function R = Resample(Expt, nresample, varargin)
%resampe an Expt


rcargs = {};
setslices = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'slices',5)
        j = j+1;
        setslices = varargin{j};
    else
        rcargs = {rcargs{:} varargin{j}};
    end
    j = j+1;
end

if Expt.Header.rc
    [rc,b] = PlotRevCorAny(Expt, rcargs{:},'nmin',1);
    nt = length(Expt.Trials);
    newtid = unidrnd(nt,nresample,nt);
    E = Expt;
    if isempty(setslices)
        slices = rc.resptime(2);
    else
        slices = setslices;
    end
    R.slices = slices;
    for j = 1:nresample
        tid = newtid(j,:);
        E.Trials = Expt.Trials(tid);
        [a,b] = PlotRevCorAny(E, rcargs{:},'slice',slices,'nmin',1);
        R.resps(j,:) = a.y(:);
        if isfield(a.sdfs,'counts')
            R.counts{j} = a.sdfs.counts;
        end
    end
    R.newtid = newtid;
    if length(slices) > 1
        nx = round(size(R.resps,2)./length(slices));
%dimensions of resps are resample, stimval, timeslice
        R.resps = reshape(R.resps,[nresample nx length(slices)]);       
    end
end