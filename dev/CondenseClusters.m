function [CC, GM] = CondenseClusters(C, varargin)
%[CC, GM] = CondenseClusters(C, varargin)
%strip GM fit objects out of cluster structure for PlotClusters/AllVPcs
%CondenseClusters(C, 'method', 2)
%            produces a 2-D Cell matrix 
%CondenseClusters(C, 'method', 3)
%            produces a 2-D structure matrix
%CondenseClusters(C, 'runtest') uses all three methods and times the
%get/set loop


method = 1;
rmf = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'method',5)
        j = j+1;
        method = varargin{j};
    elseif strncmpi(varargin{j},'rmfields',5)
        j = j+1;
        rmf = varargin{j};
    elseif strncmpi(varargin{j},'runtest',5)
        Clusters = C;
        f = figure;
        drawnow;
        if length(varargin) > j
            args = varargin(j+1:end);
        else
            args = {};
        end
        fprintf('Testing Method 1:\n');
        B = CondenseClusters(Clusters,args{:});
        setappdata(f,'Clusters',B);
        tic; B = getappdata(f,'Clusters'); B{1}{1}.ctime = now; setappdata(f,'Clusters',B); toc
        tic; B = getappdata(f,'Clusters'); B{1}{1}.ctime = now; setappdata(f,'Clusters',B); toc
        tic; B = getappdata(f,'Clusters'); B{1}{1}.ctime = now; setappdata(f,'Clusters',B); toc


        fprintf('Testing Method 2:\n');
        B = CondenseClusters(Clusters,'method',2);
        setappdata(f,'Clusters',B);
        tic; B = getappdata(f,'Clusters'); B{1,1}.ctime = now; setappdata(f,'Clusters',B); toc
        tic; B = getappdata(f,'Clusters'); B{1,1}.ctime = now; setappdata(f,'Clusters',B); toc
        tic; B = getappdata(f,'Clusters'); B{1,1}.ctime = now; setappdata(f,'Clusters',B); toc

        fprintf('Testing Method 2:\n');
        B = CondenseClusters(Clusters,'method',3);
        setappdata(f,'Clusters',B);
        tic; B = getappdata(f,'Clusters'); B(1,1).ctime = now; setappdata(f,'Clusters',B); toc
        tic; B = getappdata(f,'Clusters'); B(1,1).ctime = now; setappdata(f,'Clusters',B); toc
        tic; B = getappdata(f,'Clusters'); B(1,1).ctime = now; setappdata(f,'Clusters',B); toc
        return;
    end
    j =j+1;
end

if method == 3
    for j = 1:length(C)
        for p = 1:length(C{j})
            f = fields(C{j}{p});
            for k = 1:length(f)
                if isobject(C{j}{p}.(f{k}))
                    GM(j,p).(f{k}) = C{j}{p}.(f{k});
                else
                    CC(j,p).(f{k}) = C{j}{p}.(f{k});
                end
            end
        end
    end
elseif method == 2
    for j = 1:length(C)
        for p = 1:length(C{j})
            f = fields(C{j}{p});
            CC{j,p} = C{j}{p};
            for k = 1:length(f)
                if isobject(C{j}{p}.(f{k}))
                    GM(j,p).(f{k}) = C{j}{p}.(f{k});
                    CC{j,p} = rmfield(CC{j,p},f{k});
                end
            end
        end
    end
else %method 1. Best. Keep original structure, just remove GM objecsts
    for j = 1:length(C)
        for p = 1:length(C{j})
            f = fields(C{j}{p});
            CC{j}{p} = rmfields(C{j}{p},rmf);
            for k = 1:length(f)
                if isobject(C{j}{p}.(f{k}))
                    GM(j,p).(f{k}) = C{j}{p}.(f{k});
                    CC{j}{p} = rmfield(CC{j}{p},f{k});
                end
            end
            if isfield(CC{j}{p},'next')
                for c = 1:length(CC{j}{p}.next)
                    if ~isempty(CC{j}{p}.next{c})
                        f = fields(C{j}{p}.next{c});
                        for k = 1:length(f)
                        if isobject(C{j}{p}.next{c}.(f{k}))
                            GM(j,p).next{c}.(f{k}) = C{j}{p}.next{c}.(f{k});
                            CC{j}{p}.next{c} = rmfield(CC{j}{p}.next{c},f{k});
                        end
                        end
                    end
                end
            end
        end
    end
end
return;

allfields = {};
for j = 1:length(C)
    if ~isempty(C{j})
    f = fields(C{j});
    allfields = unique({allfields{:} f{:}});
    end
end
for j = 1:length(C)
    if ~isempty(C{j})
    f = fields(C{j});
        for k = 1:length(f)
            CC(j).(f{k}) = C{j}.(f{k});
        end
    end
end
if go == 0
    return;B = Ci
end
C = CC;

function C= CondenseClustersB(C, go, varargin)

allfields = {};
for j = 1:length(C)
    if ~isempty(C{j})
    f = fields(C{j});
    allfields = unique({allfields{:} f{:}});
    end
end
for j = 1:length(C)
    if ~isempty(C{j})
    f = fields(C{j});
        for k = 1:length(f)
            CC(j).(f{k}) = C{j}.(f{k});
        end
    end
end
if go == 0
    return;
end
C = CC;


function S = gm2struct(gm)

S.mu = gm.mu;
S.mixp = gm.PComponents;
S.Sigma = gm.Sigma;


function [C, GM] = CondenseA(C)


for j = 1:length(C)
    if iscell(C{j})
        [C{j} GM{j}] = CondenseA(C{j});
    else
        f = fields(C{j});
        for k = 1:length(f)
            if isobject(C{j}.(f{k}))
                GM{j} = C{j}.(f{k});
                C{j} = rmfield(C{j},f{k});
            end
        end
    end
end
