function CompareClusters(Ca,Cb, varargin)
% prints list of differnces between two Cluster structures
% CompareClusters(Ca,Cb)
verbose = 1;
plotmean = 1;
labels = {};
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'labels',5)
    j = j+1;
    labels = varargin{j};
 
elseif strncmpi(varargin{j},'quiet',5)
    verbose = 0;
end
j = j+1;
end

fa = fieldnames(Ca);
fb = fieldnames(Cb);
f = intersect(fa,fb);
for j = 1:length(f)
    if isstruct(Ca.(f{j}))
    elseif iscell(Ca.(f{j}))
    elseif isobject(Ca.(f{j}))
    elseif sum(size(Ca.(f{j}))) ~= sum(size(Cb.(f{j})))        
    elseif strcmp(f{j},'DprimeUsed')
    elseif sum(strcmp(f{j},{'TemplateUsed' 'mumeanUsed'}))
    elseif strcmp(f{j},'savetime')
        if verbose
        fprintf('%s %s %s\n',f{j},datestr(Ca.(f{j})(end)),datestr(Cb.(f{j})(end)));
        end
    elseif sum(Ca.(f{j}) ~= Cb.(f{j}))
        if verbose
        if ischar(Ca.(f{j}))
            fprintf('%s %s %s\n',f{j},Ca.(f{j}),Cb.(f{j}));
        else
            fprintf('%s%s::%s\n',f{j},sprintf(' %.2f',Ca.(f{j})),sprintf(' %.2f',Cb.(f{j})));
        end
        end
    end
end

if verbose == 0 %%quiet, not silent
    dropi(1) = Ca.dropi(3);
    dropi(2) = Cb.dropi(3);
    if length(labels) == 2
        s{1,1} = 'C';
        s{2,1} = labels{1};
        s{3,1} = labels{2};
        c = 1;
    elseif length(labels) ==3
        s{1,1} = labels{1};
        s{2,1} = labels{2};
        s{3,1} = labels{3};
        c = 1;        
    else
        c = 0;
    end
    s{1,1+c} = 'Events';
    s{1,2+c} = 'Cl1';
    s{1,3+c} = 'Dropi';
    s{1,4+c} = 'Mahal2D';
    s{1,5+c} = 'Fitdp';
    s{1,6+c} = 'User';
    s{1,7+c} = 'Date';
    cstrs = {'nspks' 'ncut' 'dropi' 'mahal2D' 'fitdprime' 'user' 'date'};
    s(2,(c+1):c+length(cstrs)) = ClusterStrings(Ca, cstrs);
    s(3,(c+1):c+length(cstrs)) = ClusterStrings(Cb, cstrs);
    nl = 3;
    for j = 1:min([length(Ca.next) length(Cb.next)]);
        if isfield(Ca.next{j},'space')
            nl = nl+1;
            s(nl,(c+1):c+length(cstrs)) = ClusterStrings(Ca.next{j}, cstrs);
            s{nl,1} = sprintf('Cl%d',Ca.next{j}.cluster);
        end
        if  isfield(Cb.next{j},'space')
            nl = nl+1;
            s(nl,(c+1):c+length(cstrs)) = ClusterStrings(Cb.next{j}, cstrs);
            s{nl,1} = sprintf('Cl%d',Cb.next{j}.cluster);
        end
    end
    StrTable(s);
end
    

if plotmean
GetFigure('Templates');
dy = max(Ca.MeanSpike.ms(:))./100;
hold off;
plot(Ca.MeanSpike.ms','r');
hold on;
plot(Cb.MeanSpike.ms'+dy,'g--');
plot(Ca.TemplateUsed','m');
plot(Cb.TemplateUsed'+dy,'b--');
end