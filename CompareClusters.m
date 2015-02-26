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
elseif strncmpi(varargin{j},'print',5)
    verbose = 0;
 
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
    elseif sum(Ca.(f{j})(:) ~= Cb.(f{j})(:))
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
    nf = 1;
    s{1,nf+c} = 'Events'; nf = nf+1;
    s{1,nf+c} = 'Cl1'; nf = nf+1;
    s{1,nf+c} = 'Dropi';nf = nf+1;
    s{1,nf+c} = 'Mahal2D';nf = nf+1;
    s{1,nf+c} = 'Fitdp';nf = nf+1;
    s{1,nf+c} = 'Trigger';nf = nf+1;
    s{1,nf+c} = 'Trigdt';nf = nf+1;
    s{1,nf+c} = 'recluster';nf = nf+1;
    s{1,nf+c} = 'User';nf = nf+1;
    s{1,nf+c} = 'Date';nf = nf+1;
    cstrs = {'nspks' 'ncut' 'dropi' 'mahal2D' 'fitdprime' 'Trigger' 'trigdt' 'recluster' 'user' 'date' };
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
dy = max(Ca.MeanSpike.ms(:))./50;
chspk = Ca.chspk;
hold off;
a = plot(Ca.MeanSpike.ms','r--','linewidth',2);
h(1) = a(1);
hold on;
a = plot(Cb.MeanSpike.ms'+dy,'g');
h(2) = a(1);
a = plot(Ca.TemplateUsed','m--','linewidth',2);
h(3) = a(1);
a = plot(Cb.TemplateUsed'+dy,'b');
h(4) = a(1);
legend(h, {'Spk1' 'Spk2' 'Tmpl1' 'Tmpl2'}) 
end