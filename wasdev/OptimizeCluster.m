function C = OptimizeCluster(DATA)

C = OptimizeClusterBoundary(DATA);
function C = OptimizeClusterBoundary(DATA)
    
    
    if DATA.currentcluster > 1
        preC = DATA.cluster.next{DATA.currentcluster-1};
    else
       preC = DATA.cluster;
    end
    
    
    if preC.shape == 0
        C  = OptimizeEllipse(DATA);
        C.crit = 1;
        C.sign = 0;
        return;
    end
    if preC.shape == 1 || preC.shape == 2
        [C, details]  = OptimizeLine(DATA);
        if isfield(DATA,'energy')
        C.sign = CheckSign(C,details.r,DATA.energy);
        else
            C.sign = C.sign;
        end
        return;
    end
    usegm = 1;
    
    space = preC.space;
    if space(1) == 5
        C = OptimizeVarE(DATA);
        return;
    elseif space(1) == 6  %% cut in > 2 dimensions.
        if DATA.currentcluster == 1 %should always be tru
            x = DATA.ndxy(:,1);
            y = DATA.ndxy(:,2);
        else
            x = DATA.xy{DATA.currentcluster}(:,1);
            y = DATA.xy{DATA.currentcluster}(:,2);
        end
    else
    xi = space(2);
    yi = space(3);
    end
    if space(1) == 1
        x = DATA.pcs(:,xi);
        y = DATA.pcs(:,yi);
    elseif space(1) == 2
        xi = space(3);
        yi = space(5);
        AllV = getappdata(DATA.toplevel,'AllV');
        x = squeeze(AllV(space(2),xi,:));
        y = squeeze(AllV(space(4),yi,:));
    elseif space(1) == 3
        x = DATA.TemplateScores(:,xi);
        y = DATA.TemplateScores(:,yi);
    end
    C.y = minmax(y);
    
    if usegm
        xy = cat(2,x,y);
        [a,b,c] = BestAngleGM(xy, [],[]);
        C.angle = a;
        C.crit = c.crit;
        C.sign = CheckSign(C,c.xy(:,1),DATA.energy);
        return;
    end
    a = -pi/2:pi/36:pi/2; % use this range because this is what atan returns;
    
    for j = 1:length(a)
        xy = xyrotate(x,y,a(j));
        dip(j) = HartigansDipTest(sort(xy(:,1)));
        [aa,bb] = FindDip(xy(:,1),DATA.energy(1,:));
        mydip(j) = bb.dipsize(1);
        coeff(j) = BimodalCoeff(xy(:,1),1.5);
    end
    if space(1) == 2
        [dipval,b] = max(dip);
        dipval = coeff(b);
    else
        [dipval,b] = max(coeff);
    end
    C.angle = a(b);
    xy = xyrotate(x,y,C.angle);
    

    [crit,b] = FindDip(xy(:,1),DATA.energy(1,:),'plot',DATA.tag.dips,'gmix');
    C.crit = crit(1);
    C.crit = mean(crit([5 6])); %based on GM fit
    C.sign = b.sign;
    SetFigure(DATA.tag.dips, DATA)
    hold off;
    plot(a,dip./max(dip));
    hold on;
    plot(a,mydip./max(mydip),'r');
    plot(a,coeff./max(coeff),'g');
 
    
    
  function C = OptimizeEllipse(DATA)
     c = DATA.currentcluster;
     state.cluster = c;
     if DATA.currentcluster > 1
         guess(1:4) = DATA.cluster.next{c-1}.xyr;
         guess(5) = DATA.cluster.next{c-1}.angle;
%         guess(5) = 0;  %DATA.xy{} is rotated by ellipse angle already
         if isfield(DATA.cluster.next{c-1},'aspectratio')
             state.aspectratio = DATA.cluster.next{c-1}.aspectratio;
         else
             state.aspectratio = 1;
         end
         C = DATA.cluster.next{c-1};
     else
         C = DATA.cluster;
         guess(1:4) = DATA.cluster.xyr;
         guess(5) = DATA.cluster.angle;
%         guess(5) = 0;  %DATA.xy{} is rotated by ellipse angle already
         if isfield(DATA.cluster,'aspectratio')
             state.aspectratio = DATA.cluster.aspectratio;
         else
             state.aspectratio = 1;
         end
     end
     if isfield(DATA,'xy')
     C.xy = minmax(DATA.xy{c});
     end
setappdata(DATA.toplevel,'fitparams',[]);
setappdata(DATA.toplevel,'fithists',[]);

maxiter = 100;
state.mintype = 1; %eucliean sum to weight smaller
state.mintype = 2; %use whichever is smallest
state.mintype = 3;
state.initialar = guess(4)./guess(3);
state.rx = guess(3);
state.ry = guess(4);
state.angle = guess(5);
if isfield(DATA,'xy')
    state.xy = DATA.xy(state.cluster);
else
    state.xy = C.xy;
end


guess = guess([1:3 5]);
%guess = guess([1 2]);
options = optimset('MaxFunEvals',100000,'maxiter',maxiter,'display','off');
dipguess = MinimiseEllipseDip(C, DATA);
[dpb, bfits] = MinimiseEllipseb(dipguess, DATA,state);
[dpa, afits] = MinimiseEllipseb(guess, DATA,state);
if dpb < dpa
    fprintf('Changing starting X,Y\n');
    guess = dipguess;
end
nloop = 0;
while min(afits.r) > 1 && nloop < 10 %no events in cluster
    cc = mean(DATA.xy{state.cluster});
    dvec = guess(1:2)-cc;
    guess(1:2) = guess(1:2)-dvec * min(afits.r)/4;
    [dpa, afits] = MinimiseEllipseb(guess, DATA,state);
    nloop = nloop + 1;
end
[fittedparams,fval,exitflag, output] = fminsearch(@MinimiseEllipseb,guess,options,DATA,state);
[dp, fits] = MinimiseEllipseb(fittedparams, DATA,state);
if length(guess) == 2
    C.xyr(1:2) = fittedparams(1:2);
elseif length(guess) == 4
    C.xyr = fittedparams(1:4);
    C.xyr(4) = C.xyr(3).*state.initialar;
    C.angle = fittedparams(4);
else
    C.xyr = fittedparams(1:4);
    C.angle = fittedparams(5);
end

function [SSD,dipr, details ] = EllipseDip(params, DATA, state)
%[SSD, details ] = MinimiseEllipseb(params, DATA, state)
%find ellipze that maximizes separation of 2 Gaussian fits in 1d
%does not allow change to ellipse aspectratio
cx = params(1);
cy = params(2);
if length(params) == 2
    rx = state.rx;
    ry = state.ry;
    a = state.angle;
elseif length(params) == 4
rx = params(3);
ry = params(3) * state.initialar;
a = params(4);
end
xys = xyrotate(state.xy(:,1)-cx,(state.xy(:,2)-cy) ./state.aspectratio,a);
r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;
a = MyDip(sqrt(r));
SSD = a.dipsize(2);
dipr = a.x(a.dip(2));
details.r = sqrt(r);


function guess = MinimiseEllipseDip(C, DATA)
state.initialar = C.xyr(4)./C.xyr(3);
state.aspectratio = C.aspectratio;
state.rx = C.xyr(3);
state.ry = C.xyr(4);
state.angle = C.angle;
state.cluster = DATA.currentcluster;
if isfield(DATA,'xy')
    state.xy = DATA.xy(state.cluster);
else
    state.xy = C.xy;
end
guess(1:3) = C.xyr(1:3);
guess(4) = C.angle;
xs = -guess(3)*2:guess(3)/2:guess(3)*2;
ys = -C.xyr(4)*2:C.xyr(4)/2:C.xyr(4)*2;
rs = [C.xyr(3)/2 C.xyr(3) C.xyr(3) .* 2];

for l = 1:length(rs)
    guess(3) = rs(l);
for j = 1:length(xs)
    for k = 1:length(ys)
        guess(1) = C.xyr(1)+xs(j);
        guess(2) = C.xyr(2)+ys(k);
        [dips(j,k,l) dipr(j,k,l)] = EllipseDip(guess, DATA, state);    
    end
end
end
[a,b] = min(dips(:));
[j,k,l] = ind2sub(size(dips),b);
if isfield(DATA,'xy')
    state.xy = DATA.xy(state.cluster);
else
    state.xy = C.xy;
end
        guess(1) = C.xyr(1)+xs(j);
        guess(2) = C.xyr(2)+ys(k);
        guess(3) = rs(l);
   [a,b,c] =     EllipseDip(guess, DATA, state);
        guess(3) = rs(l).*dipr(j,k);
   [a,b,c] =     EllipseDip(guess, DATA, state);
        
function [SSD, details ] = MinimiseEllipse(params, DATA, state)
%[SSD, details ] = MinimiseEllipse(params, DATA, state)
%find ellipze that maximizes separation of 2 Gaussian fits in 1d


cx = params(1);
cy = params(2);
rx = params(3);
ry = params(4);
a = params(5);
xy = DATA.xy{state.cluster};
xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./state.aspectratio,a);
r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;
%r = Rprime(r);
details.r = r;

C = DATA.cluster;
C.clst(r < 1) = 2;
C.clst(r>=1) = 1;
[dp, fits] = Fit2Gauss(C, r, DATA);
%don't just fit on traditional dprime.  If noise is far away, can
%pay to find samll sd for cluster, even by putting boundary in the
%middle of the cluster

if isempty(fits{1}) || isempty(fits{2})
    SSD = 1e10;
else
dpa = (1 - fits{1}.mean)./fits{1}.sd;
dpb = (fits{2}.mean-1)./fits{2}.sd;
details.fits = fits;
%don't just sum. Favor the smaller number, If one dp is large, don' t let improvements
%in that swamp the smaller one
if state.mintype == 1
    dp = -(sqrt(dpa)+sqrt(dpb)).^2;
else
    dp = -(dpa+dpb);
end
if dpa < 0 || dpb < 0
    SSD = 0;
else
    SSD = dp; 
end
end
fitparams = getappdata(DATA.toplevel,'fitparams');
fitparams = [fitparams; [params dp]];
setappdata(DATA.toplevel,'fitparams',fitparams);
rhist = histc(r,linspace(0, max(r)));
fithists = getappdata(DATA.toplevel,'fithists');
fithists = cat(2,fithists,rhist);
setappdata(DATA.toplevel,'fithists',fithists);

function [SSD, details ] = MinimiseEllipseb(params, DATA, state)
%[SSD, details ] = MinimiseEllipseb(params, DATA, state)
%find ellipze that maximizes separation of 2 Gaussian fits in 1d
%does not allow change to ellipse aspectratio
cx = params(1);
cy = params(2);
if length(params) == 2
    rx = state.rx;
    ry = state.ry;
    a = state.angle;
elseif length(params) == 4
rx = params(3);
ry = params(3) * state.initialar;
a = params(4);
end
xys = xyrotate(state.xy(:,1)-cx,(state.xy(:,2)-cy) ./state.aspectratio,a);
r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;
%r = Rprime(r);
details.r = r;

C = DATA.cluster;
C.clst(r < 1) = 2;
C.clst(r>=1) = 1;
[dp, fits] = Fit2Gauss(C, r, DATA);
%don't just fit on traditional dprime.  If noise is far away, can
%pay to find samll sd for cluster, even by putting boundary in the
%middle of the cluster

if isempty(fits{1}) || isempty(fits{2})
    SSD = 1e10;
else
dpa = (1 - fits{1}.mean)./fits{1}.sd;
dpb = (fits{2}.mean-1)./fits{2}.sd;
details.fits = fits;
%
%if mintype == 1 don't just sum. Favor the smaller number, If one dp is large, don' t let improvements
%in that swamp the smaller one
%if mintype == 2, just use whichever dp is the worst. Otherwise can do
%funny things to make dp for inside ellispe very good, and allow dp for
%outsiders to be much owrks
if state.mintype == 1
    dp = -(sqrt(dpa)+sqrt(dpb)).^2;
elseif state.mintype == 2
    dp = -min([dpa./sqrt(fits{1}.amp) dpb./sqrt(fits{2}.amp)]);
elseif state.mintype == 3
      if ~isfield(fits{1},'fitted') || ~isfield(fits{2},'fitted')
          dp = 1e6;
      else
          dp = fits{1}.fitted(end)./sqrt(fits{1}.amp)+fits{2}.fitted(1)./sqrt(fits{2}.amp);
      end
 else
    dp = -(dpa+dpb);
end
if dpa < 0 || dpb < 0
    if state.mintype == 3
        SSD = 1e6;
    else
        SSD = 0;
    end
else
    SSD = dp; 
end
end
if abs(imag(SSD)) > abs(SSD)/1000 %shouldn't happen
    SSD
end
fitparams = getappdata(DATA.toplevel,'fitparams');
fitparams = [fitparams; [params dp]];
setappdata(DATA.toplevel,'fitparams',fitparams);
rhist = histc(r,linspace(0, max(r)));
fithists = getappdata(DATA.toplevel,'fithists');
fithists = cat(2,fithists,rhist);
setappdata(DATA.toplevel,'fithists',fithists);

 function [C, fits] = OptimizeLine(DATA)
     c = DATA.currentcluster;
     state.cluster = c;
     if DATA.currentcluster > 1
         guess(1:4) = DATA.cluster.next{c-1}.xyr;
         guess(5) = DATA.cluster.next{c-1}.angle;
%         guess(5) = 0;  %DATA.xy{} is rotated by ellipse angle already
         if isfield(DATA.cluster.next{c-1},'aspectratio')
             state.aspectratio = DATA.cluster.next{c-1}.aspectratio;
         else
             state.aspectratio = 1;
         end
         C = DATA.cluster.next{c-1};
     else
         C = DATA.cluster;
         guess(1) = DATA.cluster.crit;
         guess(2) = DATA.cluster.angle;
         state.aspectratio = 1;
     end
if isfield(DATA,'xy')
    xy = DATA.xy{state.cluster};
else
    xy = C.xy;
end
setappdata(DATA.toplevel,'fitparams',[]);
setappdata(DATA.toplevel,'fithists',[]);
%First find best starting point - fitting routine likes local minima
state.mintype = 1; %eucliean sum to weight smaller
[a,b] = GMDip(xy,[]);
crits = [a(1) a(2) DATA.cluster.crit];
setappdata(DATA.toplevel,'fitparams',[]);
for j = 1:length(crits)
    guess(1) = crits(j);
    [dps(j), afits] = MinimiseLine(guess, DATA,state);
end
dip = MyDip(xy(:,1));
j = j+1;
guess(1) = dip.dip(2);
[dps(j), afits] = MinimiseLine(guess, DATA,state);
crit(j) = guess(1);

[a,besti] = min(dps);
guess(1) = crits(besti);

maxiter = 100;
options = optimset('MaxFunEvals',100000,'maxiter',maxiter,'display','off');
nloop = 0;
if min(afits.counts) ==0 
    [d, details] = GMDip(afits.r,[]);
    guess(1) = d(1);
end
%first optimize line position in1D
%setappdata(DATA.toplevel,'fitparams',[]);
%[fittedparams,fval,exitflag, output] = fminsearch(@MinimiseLine,guess(1),options,DATA,state);
%guess(1) = fittedparams(1);
setappdata(DATA.toplevel,'fitparams',[]);
[fittedparams,fval,exitflag, output] = fminsearch(@MinimiseLine,guess,options,DATA,state);
[dp, fits] = MinimiseLine(fittedparams, DATA,state);
C.crit = fittedparams(1);
C.angle = fittedparams(2);
C.fitdprime = fits.dpa+fits.dpb;
C.y = minmax(fits.y);


function Plot2GaussFit(params, DATA, state)

    if length(params) <= 2 %line
        C = DATA.cluster;
        C.crit = params(1);
        xy = DATA.xy{state.cluster};
        if length(params) == 2
            a = params(2);
            xys = xyrotate(xy(:,1),xy(:,2),a-C.angle);
        else
            xys = xy;
        end
        
        r = xys(:,1);
    end
    [dp, fits] = Fit2Gauss(C, r, DATA,'plothist')

function [SSD, details ] = MinimiseLine(params, DATA, state)

C = DATA.cluster;
if isfield(DATA,'xy')
    xy = DATA.xy{state.cluster};
else
    xy = C.xy;
end
if length(params) == 2
a = params(2);
xys = xyrotate(xy(:,1),xy(:,2),a-C.angle);
else
    xys = xy;
end
C.crit = params(1);
crit = C.crit;
r = xys(:,1);
details.r = r;
C.clst(r < C.crit) = 2;
C.clst(r>=C.crit) = 1;
details.counts = [sum(r < C.crit) sum(r >= C.crit)];

[dp, fits] = Fit2Gauss(C, r, DATA);
%don't just fit on traditional dprime.  If noise is far away, can
%pay to find samll sd for cluster, even by putting boundary in the
%middle of the cluster

if isempty(fits{1}) || isempty(fits{2})
    SSD = 1e10;
else
dpa = (crit - fits{1}.mean)./fits{1}.sd;
dpb = (fits{2}.mean-crit)./fits{2}.sd;
%dpa and dpb are both positive if fit means are the
%correct side of boundary. If this is not true, the value is
%uninterpretable
details.fits = fits;
details.y = xys(:,2);
details.dpa = dpa;
details.dpb = dpb;
%don't just sum. Favor the smaller number, If one dp is large, don' t let improvements
%in that swampt th esmaller one
if state.mintype == 1
    if dpa > 0 &&  dpb > 0
        dp = -(sqrt(dpa)+sqrt(dpb)).^2;
    elseif dpa > 0
        dp = -dpa/10;
    else
        dp = -dpb/10;
    end
else
    if dpa > 0 &&  dpb > 0
        dp = -(dpa+dpb);
    elseif dpa > 0
        dp = -dpa/10;
    else
        dp = -dpb/10;
    end
end
SSD = dp;
end
fitparams = getappdata(DATA.toplevel,'fitparams');
fitparams = [fitparams; [params dp]];
setappdata(DATA.toplevel,'fitparams',fitparams);


function [dp, fits, details] = Fit2Gauss(C, r, DATA, varargin)
%[dp, fits, details] = Fit2Gauss(C, r, DATA, varargin)
%negative values of dp are good separation
    plottype = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'plothist',6)
            plottype = varargin{j};
        end
        j = j+1;
    end
    fx = linspace(min(r),max(r),200);
    if isfield(C,'clst')
        id = find(C.clst == DATA.currentcluster+1);
        nid = find(C.clst ==1);
    else
    id = find(DATA.clst == DATA.currentcluster+1);
    nid = find(DATA.clst ==1);
    end
    if C.shape(1) == 0 
        r = Rprime(r);
        fx = linspace(min(r),max(r),200);
        [y,x] = hist(r,fx);
        [a,b] = min(abs(fx-1));
    elseif isfield(C,'crit')
        [a,b] = min(abs(fx-C.crit));
        [y,x] = hist(r,fx);
    else
        [a,b] = min(abs(fx-mean(fx))); %temporary
        [y,x] = hist(r,fx);
    end
    crit = x(b);
    dp = abs((mean(y(1:b)) - x(b))./std(y(1:b)));
    guess = [mean(y(1:b)) std(y(1:b)) prctile(y(1:b),95)];
    if length(id) <= 1 || length(nid) <= 1
        dp = NaN;
        fits{1}.params = [NaN NaN NaN];
        fits{1}.mean = NaN;
        fits{1}.sd = NaN;
        fits{2}.params = [NaN NaN NaN];
        fits{2}.mean = NaN;
        fits{2}.sd = NaN;
        details.fitpos = [0 0];
        return;
    end
    fits{1} = FitGauss(x(1:b),y(1:b));
    fits{2} = FitGauss(x(b:end),y(b:end));
     if isfield(fits{1},'params') && isfield(fits{2},'params')
         details.fitpos = [fits{1}.params(1) < fx(b) fits{2}.params(1) > fx(b)];
%dont need abs fit 1 ix for x < b, fit 2 is x > b. If mean 1 > mean 2, dprime is bad         
         dp = (fits{1}.params(1)-fits{2}.params(1))./sqrt((fits{1}.params(2).^2+fits{2}.params(2).^2)/2);
         if sum(details.fitpos) < 2  %one mean is wrong side of criterion. Must be bad
             dp = abs(dp);
         end
         details.diff = abs((fits{1}.params(1)-fits{2}.params(1)));
         details.sd = sqrt((fits{1}.params(2).^2+fits{2}.params(2).^2)/2);
         if strcmp(plottype,'plothist')
             GetFigure('FitGauss');
             hold off;
             bar(x,y);
             hold on;
             bar(x(1:b),y(1:b),'r')
             fya = fitGauss(fx, fits{1}.params, 'eval');
             fyb = fitGauss(fx, fits{2}.params, 'eval');
             plot(fx,fya,'r-','linewidth',2);
             plot(fx,fyb,'b-','linewidth',2);
             set(gca,'ylim',[0 max(y)]);
             title(sprintf('Dprime from fits %.2f',dp));
         end
     else
         if isempty(fits{1})
             fits{1}.params = [NaN NaN NaN];
             fits{1}.mean = NaN;
             fits{1}.sd = NaN;
         end
         if isempty(fits{2})
             fits{2}.params = [NaN NaN NaN];
             fits{2}.mean = NaN;
             fits{2}.sd = NaN;
         end
         dp = NaN;
         details.fitpos = [0 0];
         details.diff = NaN;
         details.sd = NaN;
     end

function x = Rprime(r)
    x = sqrt(r);
