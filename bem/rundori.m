function R = rundori(varargin)

R.oris = 0:pi/32:pi/4;
nruns = 100;
do = atan(0.1);
dsf = atan(0.1);
dgs = [-1:0.1:1];
dfs = [-0.5:0.05:0.5];
disps = [-0.4:0.05:0.4];
%for an sd of 0.4 a gradient of 1 is equalant do a disparity of 0.4
type = 0;
showprogress = 1;
args = {};

j = 1;
while j <= length(varargin) 
    if strncmpi(varargin{j},'nruns',3)
        j = j+1;
        nruns = varargin{j};
    elseif strncmpi(varargin{j},'ori',3)
        j = j+1;
        R.oris = varargin{j};
    elseif strncmpi(varargin{j},'dori',3) %Rf dori
        j = j+1;
        do= varargin{j};
        type = 1;
        if length(do) > 1
            type = 2;
        end
    elseif strncmpi(varargin{j},'dsf',3) %RF dSF
        j = j+1;
        dsf= varargin{j};
        type = 3;
        if length(dsf) > 1
            type = 3;
        end
    elseif strncmpi(varargin{j},'dgs',3) %disp gradient V
        j = j+1;
        dgs= varargin{j};
    elseif strncmpi(varargin{j},'dfs',3) %disp gradient H
        j = j+1;
        dfs= varargin{j};
    elseif strncmpi(varargin{j},'silent',5) %disp gradient H
        args = {args{:} 'silent'};
        showprogress = 0;
    elseif strncmpi(varargin{j},'track',5) %disp gradient H
        args = {args{:} 'track'};
        showprogress = 0;
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end

R.nruns = nruns;
R.dgs = dgs;
R.dori = do;
R.dsf = dsf;

if type == 1
for j = 1:length(R.oris)
    fprintf('Ori %.3f at %s\n',R.oris(j),datestr(now))
    R.resps(j,:) = dori('ori',R.oris(j),'nruns',nruns,'dori',do,'disp',0,'dgs',dgs,args{:});
    if showprogress
        plot(R.resps(j,:));
        title(sprintf('Ori %.3f',R.oris(j)));
        drawnow;
    end
    R.mag(j) = range(R.resps(j,:));
end
elseif type == 2
for j = 1:length(R.oris)
    fprintf('Ori %.3f at %s\n',R.oris(j),datestr(now));
    for k = 1:length(do)
    R.resps(j,k,:) = dori('ori',R.oris(j),'nruns',nruns,'dori',do(k),'disp',0,'dgs',dgs,args{:});
    if showprogress
        plot(R.resps(j,:));
        title(sprintf('Ori %.3f',R.oris(j)));
        drawnow;
    end
    end
    R.mag(j) = range(R.resps(j,:));
end

elseif type == 3
for j = 1:length(R.oris)
    fprintf('Ori %.3f at %s\n',R.oris(j),datestr(now));
    for k = 1:length(dsf)
    [R.resps(j,k,:), details] = dori('ori',R.oris(j),'nruns',nruns,'dsf',dsf(k),'disp',0,'dfs',dfs,args{:});
    if showprogress
        plot(R.resps(j,:));
        title(sprintf('Ori %.3f',R.oris(j)));
        drawnow;
    end
    end
    R.mag(j) = range(R.resps(j,:));
    R.dfs = dfs;
    R.dsf = dsf;
    R.dgs = details.dgs;
end

else
for j = 1:length(R.oris)
    fprintf('Ori %.3f at %s\n',R.oris(j),datestr(now));
    R.resps{j} = dori('ori',R.oris(j),'nruns',nruns,'dori',do,'disps',disps,'dgs',dgs,args{:});
    R.mag(j) = range(R.resps{j}(:));
    R.sresps{j} = dori('ori',R.oris(j),'dfs',dgs,'disps',disps,'nruns',nruns,'dsf',dsf,args{:});
    R.smag(j) = range(R.sresps{j}(:));
end
R.sfs = dfs;
R.disps = disps;
end



