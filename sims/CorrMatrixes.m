function [diffs, corrs, Z] = CorrMatrixes(varargin)
%[diffs, corrs, Z] = CorrMatrixes(varargin)
%make correlation matix allowing non-orthogonal kernels
%Corrmatrixes('diffs') polts mean corr as a function of diff in angle
%need to think about how this looks for C+N task.
%its only the different pool, dPD > 135 that shows reversal
%
%CorrMatrixes('manual',X) uses the correlation matrix given in X
%X = CorrMatrixes('mkmatrix',n) returns a matrix for use with 'manual'
%   n = 1: simple feature attention.
%      'nulloffset', n  moves position of second feature along diagnonal
%   n = 2: diagonal matrix (bottom up, fixed), r = cos diff
%   n = 3: like 2, but with high corr for 180 deg diff
%   n = 4: square matrix like Shadlen et al
%   n = 5: use data from Cohen et at
%   n = 6: diagonal matrix using means both conditions in Cohen et al
%   n = 7: diagonal matrix using fixation data in Cohen et al

plottype = 'matrix';
vals = 1:1140;
manualmatrix = [];
xv = 90:360:max(vals);
yv = -30:360:max(vals);
crit = abs(cos((yv(1)-xv(1)) * pi/180));
if crit < 0.01
    crit = 0.01;
end
crit = 0.9;
j = 1;
margs = {};

while j <= length(varargin)
    if strncmpi(varargin{j},'diffs',5)
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'mkmatrix',5)
        if length(varargin) > j && isnumeric(varargin{j+1})
            diffs = MakeMatrix(varargin{j+1:end});
        else
            if length(varargin) > j
                diffs = MakeMatrix(1, varargin{1+j:end});
            else
                diffs = MakeMatrix(1);
            end
        end
        return;

    elseif strncmpi(varargin{j},'manual',5)
        j = j+1;
        manualmatrix = varargin{j};
    end
    j = j+1;
end

allv = cat(2,xv,yv);
Z = [];
if isempty(manualmatrix)
for j = 1:length(allv)
for k = 1:length(allv)
    [x, y, z] = gauss2d(20,vals, 'mean', [allv(j) allv(k)]);
    a = cos((allv(j)-allv(k))*pi/180);  
    if a < crit
        z = -z;
    end
    if isempty(Z)
        Z = z;
    else
        Z = Z+z;
    end
end
end
else
    Z = repmat(manualmatrix,3,3);
    nc = size(manualmatrix,1);
    sz = floor(nc/4);
    cn = floor(nc/2);
    GetFigure('CohenNewsome');
    colors = mycolors;
    hold off;
    for p = 1:sz
        o = p-1;
    for j = 1:sz
        diffpool(p,j) = Z(o+sz+j-1,o+sz-(j-1));
        samepool(p,j) = Z(o+cn+j-1,o+cn-(j-1));
    end
    plot(samepool(p,:),'-','color',colors{p});
    hold on;
    plot(diffpool(p,:),'--','color',colors{p});
    end
    plot(mean(samepool),'k-');
    plot(mean(diffpool),'k--');
    GetFigure('Bondy')
    for j = 1:sz;
        M = TaskMatrix(nc,j);
        
    end
    return;
end

diffs = [];
corrs = [];
for j = 91:270
    for k = 91:(90+360)
        diffs = [diffs mod(j-k,360)];
        corrs = [corrs Z(j,k)];
    end
end
id = find(diffs > 180);
diffs(id) = 360-diffs(id);
hold off;
if strcmp(plottype ,'diffs')
plot(diffs,corrs,'o');
hold on;
dv = unique(diffs);
for j = 1:length(dv)
    id = find(diffs == dv(j));
    plot(dv(j),mean(corrs(id)),'ro');
end
else
imagesc(Z);
end


function M = MakeMatrix(type, varargin)


nc = 37;
nd = nc-1;
noffset = 0;
sd = 4;
task = 10;
j = 1;
plotmatrix = 1;

while j <= length(varargin)
    if strncmpi(varargin{j},'nulloff',6)
        j = j+1;
        noffset = varargin{j};
    elseif strncmpi(varargin{j},'task',4)
        j = j+1;
        task = varargin{j}; 
    end
    j = j+1;
end

M = zeros(nc);
[x,y] = meshgrid(1:nc,1:nc);
pt(1) = task; 
pt(2) = task + ceil(nc/2)-noffset;
if ismember(type,[5 6 7]) %Use Cohen+Newsome data
    corrs = [.2055 0.1593 0.1105 0.0424 0; 0.1626 0.0970 0.0684 0.1109 0]';
    corrs(:,3) = [0.2240 0.1668 0.0811 0.0731 0];
    corrs(:,4) = mean(corrs(:,1:2),2);
    T = 2-TaskMatrix(nc, task);
for j = 1:nc;
    for k = 1:nc;
        if type == 5
        c = T(j,k);
        elseif type == 6
        c = 4;
        elseif type == 7
        c = 3;
        end
        if abs(j-k) <= nd/8
            M(j,k) = corrs(1, c);
        elseif abs(j-k) > 7 * nd/8
            M(j,k) = corrs(1,c);
        elseif abs(j-k) > 6 * nd/8
            M(j,k) = corrs(2,c);
        elseif abs(j-k) > 5 * nd/8
            M(j,k) = corrs(3,c);
        elseif abs(j-k) <= nd/4
            M(j,k) = corrs(2,c);
        elseif abs(j-k) <=  3* nd/8
            M(j,k) = corrs(3,c);
        elseif abs(j-k) <=  5*nd/8
            M(j,k) = corrs(4,c);
        end
    end
end
    
elseif type == 1
    ra = abs(x + i.*y - (pt(1) +i * pt(1)));
    rb= abs(x + i.*y - (pt(1)+nc-1 + i * (pt(1)+nc-1)));
    rc= abs(x + i.*y - (pt(1) + i * (pt(1)+nc-1)));
    rd= abs(x + i.*y - (pt(1)+nc-1 + i * (pt(1))));
    r = min(cat(3,ra,rb,rc,rd),[],3);
    A = exp(-(r.^2)./(2.*sd.^2));
    M = M+A;
    r = abs(x + i.*y - (pt(2) +i * pt(2)));
    A = exp(-(r.^2)./(2.*sd.^2));
    M = M+A;
    ra = abs(x + i.*y - (pt(1) +i * pt(2)));
    rb= abs(x + i.*y - (pt(2) + i * (pt(1))));
    rc= abs(x + i.*y - (pt(2) + i * (pt(1)+nc-1)));
    rd= abs(x + i.*y - (pt(1)+nc-1 + i * (pt(2))));
    r = min(cat(3,ra,rb,rc,rd),[],3);
    A = exp(-(r.^2)./(2.*sd.^2));
    M = M-A;
    %id = sub2ind(size(M),[1 1 nc nc],[1 nc 1 nc]);
    %M(id) = 1;
    %id = sub2ind(size(M),[1 2 nc nc-1],[2 1 nc-1 nc]);
    %M(id) = 0.8;
    c = ceil(nc/2);
    %id = sub2ind(size(M),[c c c-1 c+1],[c-1 c+1 c c]);
    %M(c,c) = 1;
    %M(id) = 0.8;
    %id = sub2ind(size(M),[c 1 c nc],[1 c nc c]);
    %M(id) = -1;
elseif type == 4
M = TaskMatrix(nc, task);
else
for j = 1:nc;
    for k = 1:nc;
        M(j,k) = cos((j-k) * 2* pi/nc); 
        if abs(j-k) == 5 && type == 3
            M(j,k) = 0.5;
        end
    end
end
end

if plotmatrix
 imagesc([0:360],[0:360],M);
 set(gca,'xtick',[0 90 180 270 360],'ytick',[0 90 180 270 360]);
end

function M = TaskMatrix(nc, task)
M = zeros(nc);
for j = 1:nc;
    for k = 1:nc;
        if abs(j-task) <nc/4 && abs(k-task) < nc/4
            M(j,k) = 1;
        elseif abs(j-task-nc/2) < nc/4 && abs(k-task-nc/2) < nc/4;
            M(j,k) = 1;
        elseif abs(j-task-nc) < nc/4 && abs(k-task-nc) < nc/4;
            M(j,k) = 1;
        end
    end
end
