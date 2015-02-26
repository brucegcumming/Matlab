function [dip, details] = GMDip(xy, energy, varargin)
% [dip, details] = GMDip(xy, energy, varargin)
%fits a mixture of two Gaussians in one dimension to xy(:,1)
%
%returns  dip - three estimates of dividing line ()
%details.mahal  vector of mahal distances for different starting points
%details.G

details.type = 1;
plotdips = NaN;  %no plo tby default
niter = 100;
mydip = 0;
crit = [];
idlist = [];
dip = NaN;
usecritonly = 0;
label = '';

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'crit',4)
        if strncmpi(varargin{j},'critonly',8)
            usecritonly = 1;
        end
        j = j+1;
        crit = varargin{j};
    elseif strncmpi(varargin{j},'idlist',6)
        j = j+1;
        idlist = varargin{j};
    elseif strncmpi(varargin{j},'label',4)
        j = j+1;
        label = varargin{j};
    elseif strncmpi(varargin{j},'noplot',4)
        plotdips = NaN;
    elseif strncmpi(varargin{j},'plot',4)
        plotdips = 1;
    end
    j = j+1;
end

if length(xy) == 0
    return;
end
if diff(size(xy)) > 1
xy = xy';
end
v = sort(xy(:,1));

sigma = std(v);
if mydip
    sm = sigma/10;
    r = range(xy(:,1));
    x = min(xy(:,1)) - r/10:r/100:max(xy(:,1))+r/10;
    [y,x] = smhist(v,'sd',sm,'xval',x);
    peaks = find(diff(sign(diff(y))) < 0);
    while length(peaks) < 2
        sm = sm ./ 2;
        [y,x] = smhist(v,'sd',sm,'xval',x);
        peaks = find(diff(sign(diff(y))) < 0);
    end
    while length(peaks) > 2
        sm = sm .* 1.5;
        [y,x] = smhist(v,'sd',sm);
        peaks = find(diff(sign(diff(y))) < 0);
    end
    while length(peaks) < 2
        sm = sm * 0.9;
        [y,x] = smhist(v,'sd',sm);
        peaks = find(diff(sign(diff(y))) < 0);
    end
    peaks = peaks([1 end]);
    dip = x(peaks);
    rawx = x;
    rawy = y;
    
    a = find(y > max(y)/40);
    id = a(1):a(end);
    dpeaks = find(abs(diff(sign(diff(y(id),2)))) > 0);
    while length(dpeaks) > 4
        sm = sm  .* 2;
        y = smhist(v,'sd',sm,'xval',x);
        dpeaks = find(abs(diff(sign(diff(y(id),2)))) > 0);
    end
    signpeaks = diff(sign(diff(y(id),2)));
    signpeaks = signpeaks(dpeaks);

    
    lastbig = 0;
    mx = 0.8;
    j = 0;
    while length(dpeaks) ~= 3 && length(dpeaks) ~= 4
        j = j+1;
        if length(dpeaks) > 3
            sm = sm ./ mx;
            big = 1;
        else
            sm = sm.*mx; %less smoothing
            big = 0;
        end
        sms(j) = sm;
        if lastbig ~=big;
                mx = sqrt(mx);
        end
            lastbig = big;
        y = smhist(v,'sd',sm,'xval',x);
        dpeaks = find(abs(diff(sign(diff(y(id),2)))) > 0);
        signpeaks = diff(sign(diff(y(id),2)));
        signpeaks = signpeaks(dpeaks);
        np(j) = length(dpeaks);
%        fprintf('sm %.3f, %d peaks\n',sm,length(dpeaks));
    end
    fy = smhist(v,'sd',sm/2,'xval',x);

    dy = diff(fy(id));
    [a,b] = max(y);  %1 or 3
    ns = sum(x(id(dpeaks)) < x(b));
    if length(dpeaks) == 3
        if signpeaks(1) > 0
        [a,b] = min(abs(dy(dpeaks(1):dpeaks(2))));
        dip(4) = x(id(dpeaks(1)+b-1));
        else
        [a,b] = min(abs(dy(dpeaks(2):dpeaks(3))));
        dip(4) = x(id(dpeaks(2)+b-1));
        end
    elseif ns > 2
        [a,b] = min(abs(dy(dpeaks(2):dpeaks(3))));
        dip(4) = x(id(dpeaks(2)+b-1));
    else
        [a,b] = min(abs(dy(dpeaks(2):dpeaks(3))));
        dip(4) = x(id(dpeaks(2)+b-1));
    end
end

if isobject(energy) %evaluate a fit
    G = energy;
    x = linspace(min(v), max(v),500);
    y = pdf(G, x');
    details.gxy(:,1) = x;
    y = pdf('norm',x,G.mu(1), sqrt(G.Sigma(1))) .* G.PComponents(1);
    details.gxy(:,2) = y;
    z = pdf('norm',x,G.mu(2), sqrt(G.Sigma(2))) .* G.PComponents(2);
    details.gxy(:,3) = z;
        gid = find(x > min(G.mu) & x < max(G.mu));
    if isempty(gid)
        dip(2) = mean(G.mu);
    else
        [aa, peak] = min(z(gid)+y(gid));
        dip(2) = x(peak+gid(1)-1);
    end
    peak = find(diff(sign(z-y)) ~= 0);
    [aa,bb] = max(z(peak)+y(peak));
    details.gxy(:,3) = z;
    if isempty(bb)
        dip(3) = mean(G.mu);
    else
    dip(3) = x(peak(bb));
    end
    dip(1) = mean(dip(2:3));
    details.dipsize = 0;
    details.cdipsize = 0;
    details.G = G;
    details.gmdprime = gmdprime(G);
    details.dip = dip;
    return;
end

if usecritonly == 0
    try
    Gs{1} = gmdistribution.fit(xy(:,1),2,'Options',statset('MaxIter',niter)); %2 Gaussians, 1 dimension
    d(1) = gmdprime(Gs{1});
    details.dipres(1) = 1;
    catch
        fprintf('GMDip 1 %s fit failed %s)\n',label,lasterr);
        Gs{1}.mu = [0 0];
        Gs{1}.Converged = -1;
        d(1) = 0;
        details.dipres(1) = 0;

    end
    S.mu(1,1) = mean(v)+sigma;
    S.mu(2,1) = mean(v)-sigma;
    S.Sigma(:,:,1) = sigma.^2./2;
    S.Sigma(:,:,2) = sigma.^2/2;
    S.PComponents = [0.5 0.5];
    try 
        Gs{2} = gmdistribution.fit(xy(:,1),2,'Options',statset('MaxIter',niter),'Start',S);
        d(2) = gmdprime(Gs{2});
        details.dipres(2) = 1;
    catch
        fprintf('GMDip 2 fit (%s%.3f, %.3f sd %.3f) failed-%s\n',label,S.mu(1,1),S.mu(2,1),sigma,lasterr);
        Gs{2} = S;
        Gs{2}.Converged = -1;
        d(2) = 0;
        details.dipres(2) = 0;
    end
    nf = 2;
else
    nf = 0;
end


    for j = 1:length(crit)
        id = find(v > crit);
        nid = find(v <= crit);
        S.mu(1,1) = mean(v(id));
        S.mu(2,1) = mean(v(nid));
        S.Sigma(:,:,1) = var(v(id));
        S.Sigma(:,:,2) = var(v(nid));
        S.PComponents = [length(id) length(nid)]./length(v);
    try 
        Gs{nf+j} = gmdistribution.fit(xy(:,1),2,'Options',statset('MaxIter',niter),'Start',S);
        d(nf+j) = gmdprime(Gs{nf+j});
        details.dipres(nf+j) = 1;
    catch
        fprintf('GMDip 2 fit (%s%.3f, %.3f sd %.3f) failed %s\n',label,S.mu(1,1),S.mu(2,1),sigma,lasterr);
        Gs{nf+j} = S;
        Gs{nf+j}.Converged = -1;
        d(nf+j) = 0;
        details.dipres(nf+j) = 0;
    end
    j = length(d);
        
       if length(idlist) == size(xy,1) %give a list of classifications. Start with fit to their distributions
        j = j+1;
        try
            id = find(idlist == 2);
            nid = find(idlist==1);
        S.mu(1,1) = mean(v(id));
        S.mu(2,1) = mean(v(nid));
        S.Sigma(:,:,1) = var(v(id));
        S.Sigma(:,:,2) = var(v(nid));
        S.PComponents = [length(id) length(nid)]./length(v);
            Gs{j} = gmdistribution.fit(xy(:,1),nd,'Options',statset('MaxIter',1000),'Start',S);
            d(j) = gmdprime(Gs{j});
        catch
            fprintf('GM Fit %s(Start with Classification) fail%s\n',label,lasterr);
            Gs{j} = S;
            Gs{j}.Converged = -1;
            Gs{j}.NlogL = NaN;
            d(j) = 0;
        end
    end

    end
    
    if mydip
    S.mu(1,1) = x(peaks(1));
    S.mu(2,1) = x(peaks(2));
    S.PComponents(1) = y(peaks(1))./sum(y(peaks));
    if S.PComponents(1) > 0.9
        S.PComponents(1) = 0.8;
        id = find(y(peaks(1):peaks(2)) > max(y)./3);
        [a,b] = min(y(id+peaks(1)-1));
        S.mu(2,1) = x(id(b)+peaks(1)-1);
    end
    if S.PComponents(1) < 0.1
        S.PComponents(1) = 0.2;
        id = find(y(peaks(1):peaks(2)) > max(y)./3);
        [a,b] = min(y(id+peaks(1)-1));
        S.mu(1,1) = x(id(b)+peaks(1)-1);
    end
    S.PComponents(2) = 1 - S.PComponents(1);
    S.Sigma(:,:,1) = sigma.^2./2;
    S.Sigma(:,:,2) = sigma^2./2;
    try
        Gs{3} = gmdistribution.fit(xy(:,1),2,'Options',statset('MaxIter',niter),'Start',S);
        d(3) = gmdprime(Gs{3});
        details.dipres(3) = 1;
    catch
        fprintf('GMDip 3 fit %sfailed %s\n',label,lasterr);
        Gs{3} = S;
        Gs{3}.Converged = -1;
        d(3) = 0;
        details.dipres(3) = 0;
    end
    end
    
    for j = 1:length(Gs)
        details.converged(j) = Gs{j}.Converged;
    end
if usecritonly == 0
    if details.converged(2) == 0 && plotdips == 0
        plotdips = 1;
    end
end
    
    id = find(details.converged >= 0);
    if isempty(id)
        [D,b] = max(d);
    else
        [D,b] = max(d(id));
        b = id(b);
    end
    G = Gs{b};
    details.G = Gs;
    details.mahal = d;
    details.gmdprime = D;
    details.dipsize = 0;
    details.cdipsize = 0;
    details.best = b;

    x = linspace(min(v), max(v),500);
    if G.Converged >= 0
    y = pdf(G, x');
    details.gxy(:,1) = x;
    y = pdf('norm',x,G.mu(1), sqrt(G.Sigma(1))) .* G.PComponents(1);
    details.gxy(:,2) = y;
    z = pdf('norm',x,G.mu(2), sqrt(G.Sigma(2))) .* G.PComponents(2);
    else
        y = zeros(size(x));
        z = zeros(size(x));
    end

    if plotdips > 0
        if ~exist('rawx')
            r = range(xy(:,1));
            rawx = min(xy(:,1)) - r/10:r/100:max(xy(:,1))+r/10;
            sm = sigma/10;
            rawy = smhist(v,'sd',sm,'xval',rawx);
        end
    
        cf = gcf;
    GetFigure('Dips');
    hold off;
    [a,b] = hist(v,500);
    bar(b,a,1);
    hold on;
    scale = trapz(b,a);
    plot(x,(y+z).*scale,'r');
    plot(x,y.*scale,'g');
    plot(x,z.*scale,'g');
    sy = pdf('norm',x,S.mu(1), sqrt(S.Sigma(1))) .* S.PComponents(1);
    details.gxy(:,2) = y;
    sz = pdf('norm',x,S.mu(2), sqrt(S.Sigma(2))) .* S.PComponents(2);
    plot(x,sy.*scale,'m');
    plot(x,sz.*scale,'m');
    plot(rawx,rawy.*scale./trapz(rawx,rawy),'k');
    end


    gid = find(x > min(G.mu) & x < max(G.mu));
    
    if isempty(gid)
        dip(2) = mean(G.mu);
        peaks = max(z+y);
    else
        peaks = z(gid([1 end])) + y(gid([1 end]));
        [aa, peak] = min(z(gid)+y(gid));
        dip(2) = x(peak+gid(1)-1);
        dipsize(2) = aa./min(peaks);
    end

    peak = find(diff(sign(z-y)) ~= 0);
    [aa,bb] = max(z(peak)+y(peak));
    details.gxy(:,3) = z;
    if isempty(bb)
        dip(3) = mean(G.mu);
    else
        dip(3) = x(peak(bb));
        dipsize(3) = aa./min(peaks);
    end
    dip(1) = mean(dip(2:3));
    [aa,bb] = min(abs(dip(1)-x));
    dipsize(1) = (z(bb)+y(bb))./min(peaks);
    
    if length(energy) > 1
    e(1) = mean(energy(find(xy(:,1) > dip(1))));
    e(2) = mean(energy(find(xy(:,1) < dip(1))));
    else
        e = [0 0 ];
    end
    if mean(xy(:,1)) > std(xy(:,1))
        details.sign = 1;
    elseif mean(xy(:,1)) < -std(xy(:,1))
        details.sign = -1;
    elseif e(1) > e(2)
        details.sign = 1;
    else
        details.sign = -1;
    end
    if plotdips > 0
    plot([dip(3) dip(3)],get(gca,'ylim'),'m-');
    plot([dip(2) dip(2)],get(gca,'ylim'),'m--');
    if length(dip) > 3
    plot([dip(4) dip(4)],get(gca,'ylim'),'k');
    end
    figure(cf);
    end
    details.dipsize = dipsize;
details.dip = dip;

    
    
function d = gmdprime(G, varargin)
%calcualte drpime between two Gaussians in gmdistribution fit        
    distance = mahal(G,G.mu);
    d = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));