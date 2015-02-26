function res = RunPlotClusters(name)

res = {};
if iscellstr(name)
    for j = 1:length(name)
        x = RunPlotClusters(name{j});
        res = {res{:} x{:}};
    end
    return;
elseif iscell(name)
    CheckResults(name);
    return;
end
d = mydir(name);
for j = 1:length(d)
    PlotClusters('close');
    if d(j).isdir
        name = dir2name(d(j).name,'smrmat');
    else
        name = d(j).name;
    end
    if exist(name)
        fprintf('Running %s\n',name);
        try
            res{j} = PlotClusters(name,'load','checkrateseq');
        catch ME
            res{j}.errstate = ME;
            res{j}.filename = name;
        end
    end
end


function name = dir2name(path, type)


   [a,b] = fileparts(path);
   monkey = GetMonkeyName(path);
   name = [path '/' monkey b '.mat'];

   
function CheckResults(R)

for j = 1:length(R)
    if isfield(R{j},'datadir')
        name = [R{j}.datadir '/ExptList.mat'];
        if exist(name)
            X = load(name);
            nex = length(X.ExptList.expnames);
            [n, exlist] = Counts(X.ExptList.expnames);
        else
            X = ListExpts(R{j}.datadir);
            nex = length(X);
            [n, exlist] = Counts(X);
        end
        s = [];
        for k = 1:length(exlist)
            s = [s exlist{k} '(' num2str(n(k)) ')'];
        end
        if isfield(R{j},'exptid')
            fprintf('%s: %d cells, %d/%d Expts %s\n',R{j}.datadir,length(R{j}.cellid),length(R{j}.exptid),nex,s);
        elseif isfield(R{j},'errs')
            for k = 1:length(R{j}.errs)
                cprintf('red','%s:%d Expts %s %s\n',R{j}.datadir,nex,R{j}.errs{k},s);
            end
        else
            cprintf('blue','%s: %d cells, %d/%d Expts\n',R{j}.datadir,length(R{j}.cellid),length(R{j}.exptid),nex);
        end
    end
end