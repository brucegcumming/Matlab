function result = TestCluster(C, varargin)
%Test operatations on Cluster structs, 
%Default looks at KDE
%
plotxy = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plotxy',6)
        plotxy = 1;
    end
    j=j+1;
end


    xy = xyrotate(C.xy(:,1),C.xy(:,2),C.angle);
    if plotxy
        PlotND(C.xy,[],'idlist',C.clst);
        return; 
    end
        v = sort(xy(:,1));
    x = linspace(v(1), v(end), 200);
    dv = diff(v);
    prcs = 90:0.1:100;
    sm = linspace(prctile(dv,90),max(dv));
    for j = 1:length(sm)
        y = smhist(v,'sd',sm(j),'xval',x);
        result.minima{j} = find(diff(sign(diff(y))) > 0);
        z(j,:) = y./sum(y);
    end
    maxz = max(z');
    for j = 1:50
        th = maxz .* (j-1)./50;
        for k = 1:size(z,1)
            id = find(z(k,:) > th(j));
            dsd(j,k) = std(id);
        end
    end
    hold off;
    imagesc(x,sm,z);
    hold on;
    for j = 1:length(sm)
        if ~isempty(result.minima{j})
            plot(x(result.minima{j}),sm(j),'r.')
        end
    end
    result.z = z;
    result.sm = sm;