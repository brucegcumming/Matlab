function [gain, pos, stds] = EyeCal(Expt, varargin)
%[gain, pos, stds] = EyeCal(Expt, varargin)
%calibration of eye data from mat files made by spike 2
%EyeCal(Expt, 'summarize')
%columns are h,v,h,v,fx,fy
%rows are R, L. gains(1,1) = RH gain. gains(1,2) = gain of RV for H movement of fixpos



MEANS = 1;
VERT = 2;
MEDIANS = 3;
XYPLOT = 4;
PLOTSD = 5;
TIMEPLOT=6;

stattype = MEDIANS;
plottype = XYPLOT;
rmoffset = 0;
calcnoise =0;
mkcal = 0;
gain = [];
idlist = [];
summarize = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'allexpts',5)
        for k = 1:length(Expt.Header.ExptStart)
            if k > length(Expt.Header.ExptEnd)
                Expt.Header.ExptEnd(k) = length(Expt.Trials);
            end
            fprintf('Expt %d: %d-%d (%.1f - %.1f)\n',k,Expt.Header.ExptStart(k),Expt.Header.ExptEnd(k),...
                Expt.Trials(Expt.Header.ExptStart(k)).Start(1)/10000,Expt.Trials(Expt.Header.ExptEnd(k)).End(end)/10000);
            EyeCal(Expt,varargin{1:j-1},'expt',k);
        end
    elseif strncmpi(varargin{j},'means',4)
        stattype = MEANS;
    elseif strncmpi(varargin{j},'mkcal',4)
        mkcal =1;
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        idlist = varargin{j};
        id = find(ismember([Expt.Trials.id],idlist));
        Expt.Trials = Expt.Trials(id);
    elseif strncmpi(varargin{j},'expt',3)
        j = j+1;
        exlist = varargin{j};
        tid = [];
        for k = 1:length(exlist)
            tid = [tid Expt.Header.ExptStart(exlist(k)):Expt.Header.ExptEnd(exlist(k))];
        end
        Expt.Trials = Expt.Trials(tid);
    elseif strncmpi(varargin{j},'noise',5)
        calcnoise = 1;
    elseif strncmpi(varargin{j},'plotmean',6)
        plottype = MEANS;
    elseif strncmpi(varargin{j},'plotsd',4)
        plottype = PLOTSD;
    elseif strncmpi(varargin{j},'summary',6)
        summarize = 1;
    elseif strncmpi(varargin{j},'zoff',4)
        rmoffset = 1;
    end
    j = j+1;
end

if ~isfield(Expt.Header,'CRrates')
    Expt.Header.CRrates(1) = Expt.Header.CRsamplerate;
elseif Expt.Header.CRrates(1) < 0
    Expt.Header.CRrates(1) = 1./597;
end


if length([Expt.Trials.ftime]) ~= length([Expt.Trials.Start])
    for j = 1:length(Expt.Trials)
        if isempty(Expt.Trials(j).Start) || isempty(Expt.Trials(j).ftime)
            good(j) = 0;
        elseif Expt.Trials(j).Start < Expt.Trials(j).ftime;
            good(j) = 0;
        else
            good(j) = 1;
        end
    end
    Expt.Trials = Expt.Trials(find(good));
end
pre = mean([Expt.Trials.Start] - [Expt.Trials.ftime]);
if pre < 0
    fprintf('ERROR: Some trials start before sampling');
    pre = 0;
end
ps = round((pre./10000) ./ Expt.Header.CRrates(1));
if ps < 1
    ps = 1;
end
skip = 0;
threeeyes = 0;
addsoft = 1;
so = [0 0 0 0];

if isfield(Expt,'Comments') && ~isempty(Expt.Comments.times)
    delay = 100000; %10 sec
    if iscell(Expt.Comments.times)
        cmtimes = cat(2,Expt.Comments.times{:});
    else
        cmtimes = Expt.Comments.times;
    end
    id = find(cmtimes > Expt.Trials(1).Start(1)-delay & cmtimes < Expt.Trials(1).End(end)+delay);
    for j = 1:length(id)
        if ~strncmp(Expt.Comments.text{id(j)},'cm=rf',5) && ~strncmp(Expt.Comments.text{id(j)},'cm=noback',8) 
            fprintf('%.2f: %s\n',Expt.Comments.times(id(j)),Expt.Comments.text{id(j)});
        end
    end
end

nsd = 0;
sdlen = 100;
for j = 1:length(Expt.Trials)
    if isfield(Expt.Trials, 'softoff') && addsoft == 1;
        so = Expt.Trials(j).softoff .* 2;
    end
    
    fx(j) = Expt.Trials(j).fx;
    fy(j) = Expt.Trials(j).fy;
    if isfield(Expt.Trials,'EyeData')
            lh(j) = mean(Expt.Trials(j).EyeData(:,1))+ so(3);
            rh(j) = mean(Expt.Trials(j).EyeData(:,2))+ so(1);
            lv(j) = mean(Expt.Trials(j).EyeData(:,3))+ so(4);
            rv(j) = mean(Expt.Trials(j).EyeData(:,4))+ so(2);
            stds(j,:) = std(Expt.Trials(j).EyeData);
            rlen = size(Expt.Trials(j).EyeData,1);
            Expt.Trials(j).good = 1;
    elseif ~isempty(Expt.Trials(j).Eyevals)
        len = round(((Expt.Trials(j).End - Expt.Trials(j).Start)./10000)./Expt.Header.CRrates(1));
        rlen = min([length(Expt.Trials(j).Eyevals.lh) length(Expt.Trials(j).Eyevals.rh)...
            length(Expt.Trials(j).Eyevals.lv) length(Expt.Trials(j).Eyevals.rv)]);
        lens(j) = rlen;
        if (rlen >= len)
            sid = intersect(find(Expt.Trials(j).Eyevals.lh > -20),[skip+ps:ps+len]);
            lh(j) = mean(Expt.Trials(j).Eyevals.lh(sid))+ so(3);
            stds(j,1) = std(Expt.Trials(j).Eyevals.lh(sid));
            sid = intersect(find(Expt.Trials(j).Eyevals.lv > -20),[skip+ps:ps+len]);
            lv(j) = mean(Expt.Trials(j).Eyevals.lv(sid)) + so(4);
            stds(j,3) = std(Expt.Trials(j).Eyevals.lv(sid));
            rh(j) = mean(Expt.Trials(j).Eyevals.rh(skip+ps:ps+len)) + so(1);
            rv(j) = mean(Expt.Trials(j).Eyevals.rv(skip+ps:ps+len)) + so(2);
        stds(j,2) = std(Expt.Trials(j).Eyevals.rh(skip+ps:ps+len));
        stds(j,4) = std(Expt.Trials(j).Eyevals.rv(skip+ps:ps+len));
        ids(j) = Expt.Trials(j).id;
        if isfield(Expt.Trials(j).Eyevals,'xh')
            threeeyes = 1;
            xh(j) = mean(Expt.Trials(j).Eyevals.xh(skip+ps:ps+len));
            xv(j) = mean(Expt.Trials(j).Eyevals.xv(skip+ps:ps+len));
        end
        else
        fy(j) = NaN;
        fx(j) = NaN;
        lens(j) = 0;
        end
    else 
        lens(j) = 0;
        fy(j) = NaN;
        fx(j) = NaN;
    end
    nsd = nsd+ceil((rlen-sdlen)./10);
end    


if calcnoise
    tic;
    sds(nsd,1:6) = 0;
    nsd = 1;
    for j = 1:length(Expt.Trials)
        for k = 1:10:lens(j)-sdlen
            sds(nsd,1) = std(Expt.Trials(j).Eyevals.lh(k:k+sdlen));
            sds(nsd,2) = std(Expt.Trials(j).Eyevals.rh(k:k+sdlen));
            sds(nsd,3) = std(Expt.Trials(j).Eyevals.lv(k:k+sdlen));
            sds(nsd,4) = std(Expt.Trials(j).Eyevals.rv(k:k+sdlen));
            if threeeyes
                sds(nsd,5) = std(Expt.Trials(j).Eyevals.xh(k:k+sdlen));
                sds(nsd,6) = std(Expt.Trials(j).Eyevals.xv(k:k+sdlen));
            end
            nsd = nsd+1;
        end
    end
    toc
end
    
    
np = min([length(fx) length(fy) length(lh) length(rh)]);
good = [Expt.Trials.good] == 1; %exclude badfix trials
if np > length(good)
    np = length(good);
end
good = good(1:np);
fx = fx(1:np);
fy = fy(1:np);
lh = lh(1:np);
rh = rh(1:np);
lv = lv(1:np);
rv = rv(1:np);
pos(1,:) = lh(1:np);
pos(2,:) = rh(1:np);
pos(3,:) = lv(1:np);
pos(4,:) = rv(1:np);
pos(5,:) = fx;
pos(6,:) = fy;
aid = find(fx < 0 & ~isnan(rh+rv+lh+lv) & good);
bid = find(fx > 0 & ~isnan(rh+rv+lh+lv) & good);
fxs(1) = mean(fx(aid));
fxs(2) = mean(fx(bid));
cid = find(fy < 0  & ~isnan(rh+rv+lh+lv) & good);
did = find(fy > 0  & ~isnan(rh+rv+lh+lv) & good);
zid = find(fy == 0 & fx == 0 & good);

if stattype == MEANS
rhm(1) = mean(rh(aid));
rhm(2) = mean(rh(bid));
lhm(1) = mean(lh(aid));
lhm(2) = mean(lh(bid));
rhm(3) = mean(rh(cid));
rhm(4) = mean(rh(did));
lhm(3) = mean(lh(cid));
lhm(4) = mean(lh(did));
rvm(1) = mean(rv(cid));
rvm(2) = mean(rv(did));
lvm(1) = mean(lv(cid));
lvm(2) = mean(lv(did));
rvm(3) = mean(rv(aid));
rvm(4) = mean(rv(bid));
lvm(3) = mean(lv(aid));
lvm(4) = mean(lv(bid));
rvm(5) = mean(rv(zid));
rhm(5) = mean(rh(zid));
lvm(5) = mean(lv(zid));
lhm(5) = mean(lh(zid));
else
rhm(1) = median(rh(aid));
rhm(2) = median(rh(bid));
rhm(3) = median(rh(cid));
rhm(4) = median(rh(did));
lhm(1) = median(lh(aid));
lhm(2) = median(lh(bid));
lhm(3) = median(lh(cid));
lhm(4) = median(lh(did));
rvm(1) = median(rv(aid));
rvm(2) = median(rv(bid));
lvm(1) = median(lv(aid));
lvm(2) = median(lv(bid));
rvm(3) = median(rv(cid));
rvm(4) = median(rv(did));
lvm(3) = median(lv(cid));
lvm(4) = median(lv(did));
rvm(5) = median(rv(zid));
rhm(5) = median(rh(zid));
lvm(5) = median(lv(zid));
lhm(5) = median(lh(zid));
end
rhsd(1) = (std(rh(aid))+std(rh(bid)))/2;
rhsd(2) = mean(stds([aid bid],2));
lhsd(1) = (std(lh(aid))+std(lh(bid)))/2;
lhsd(2) = mean(stds([aid bid],1));
rvsd(1) = (std(rh(cid))+std(rh(did)))/2;
rvsd(2) = mean(stds([cid did],4));
lvsd(1) = (std(lh(cid))+std(lh(did)))/2;
lvsd(2) = mean(stds([cid did],3));

if threeeyes
    xhm(1) = median(xh(aid));
    xhm(2) = median(xh(bid));
    xhm(3) = median(xh(cid));
    xhm(4) = median(xh(did));
    xvm(1) = median(xv(aid));
    xvm(2) = median(xv(bid));
    xvm(3) = median(xv(cid));
    xvm(4) = median(xv(did));
end
fys(1) = mean(fy(cid));
fys(2) = mean(fy(did));
gain(1,1) = diff(rhm([1 2]))./diff(fxs);
gain(1,2) = diff(rvm([1 2]))./diff(fxs); %crosstalk
gain(2,1) = diff(lhm([1 2]))./diff(fxs);
gain(2,2) = diff(lvm([1 2]))./diff(fxs);
gain(1,3) = diff(rhm([3 4]))./diff(fys);
gain(1,4) = diff(rvm([3 4]))./diff(fys);
gain(2,3) = diff(lhm([3 4]))./diff(fys);
gain(2,4) = diff(lvm([3 4]))./diff(fys);
if threeeyes
gain(3,1) = diff(xhm([1 2]))./diff(fxs);
gain(3,4) = diff(xvm([3 4]))./diff(fys);
gain(3,2) = diff(xvm([1 2]))./diff(fxs);
gain(3,3) = diff(xhm([3 4]))./diff(fys);
fprintf('gain RH %.3f RV %.3f LH %.3f LV %.3f XH %.3f XV %.3f\n',gain(1,1),gain(1,4),gain(2,1),gain(2,4),gain(3,1),gain(3,4));
else
fprintf('gain RH %.3f RV %.3f LH %.3f LV %.3f\n',gain(1,1),gain(1,4),gain(2,1),gain(2,4));
fprintf('xtalk RH %.3f RV %.3f LH %.3f LV %.3f\n',gain(1,2),gain(1,3),gain(2,2),gain(2,3));

rhsd = rhsd./abs(gain(1,1)); %in real units
lhsd = lhsd./abs(gain(2,1)); %in real units
rvsd = rvsd./abs(gain(1,4)); %in real units
lvsd = lvsd./abs(gain(2,4)); %in real units
fprintf('STDs RH %.3f(%.3f) LH %.3f(%.3f) RV %.3f(%.3f) LV %.3f(%.3f)\n',...
    rhsd(1),rhsd(2),lhsd(1),lhsd(2),rvsd(1),rvsd(2),lvsd(1),lvsd(2));
end

allfy = unique(fy);
allfx = unique(fx);
n = 0;
for j = 1:length(allfx);
    for k = 1:length(allfy);
        id = find(fy == allfy(k) & fx == allfx(j));
        if length(id)
            n = n+1;
            emdat(n).fy = allfy(k);
            emdat(n).fx = allfx(j);
            emdat(n).trials = id;
            emdat(n).eyepos = [mean(rh(id)) mean(lh(id)) mean(rv(id)) mean(lv(id))];
        end
    end
end

        

if plottype == TIMEPLOT
    plot(rh,'.');
    hold on;
    if threeeyes
    plot(xh,'r.');
    else
        
    plot(lh,'r.');
    plot(lh-rh,'m.');
    end
elseif plottype == MEANS
    plot(rhm,rvm,'ro');
    hold on;
    plot(lhm,lvm,'go');
    for j = 1:length(emdat)
        plot(emdat(j).eyepos(1),emdat(j).eyepos(3),'ro','markerfacecolor','r');
        plot(emdat(j).eyepos(2),emdat(j).eyepos(4),'go','markerfacecolor','g');
    end
    
elseif plottype == XYPLOT
    allid = [aid bid cid did zid];
    tids = [Expt.Trials(allid).id];
    if rmoffset
    lhd = mean(lh(allid)-fx(allid));
    lvd = mean(lv(allid)-fy(allid));
    rhd = mean(rh(allid)-fx(allid));
    rvd = mean(rv(allid)-fy(allid));
    xvd = median(xv(allid)-fy(allid));
    xhd = median(xh(allid)-fx(allid));
    scales = [1 1 1 1];
    scales(1) = gain(2,1); %lh
    scales(2) = gain(1,1); %rh
    scales(3) = gain(2,4); %lv
    scales(4) = gain(1,4); %rv
    scales(5) = gain(3,1);
    scales(6) = gain(3,4);
    else
        lhd = 0;
        lvd = 0;
        rvd = 0;
        rhd = 0;
        xhd = 0;
        xvd = 0;
        scales = [1 1 1 1 1 1];
    end
    for j = 1:length(allid)
    hl = plot((lh(allid(j))-lhd)./scales(1),(lv(allid(j))-lvd)./scales(3),'go','buttondownfcn',{@HitPoint, tids(j)});
    hold on;
    end
    hr = plot((rh([aid bid cid did zid])-rhd)./scales(2),(rv([aid bid cid did zid])-rvd)./scales(4),'ro');
    if threeeyes
    hx = plot((xh([aid bid cid did zid])-xhd)./scales(5),(xv([aid bid cid did zid])-xvd)./scales(6),'mo');
    end
    legend([hl hr],'L','R');
    plot((mean(lh(aid))-lhd)./scales(1),(mean(lv(aid))-lvd)./scales(3),'go','markerfacecolor','g');
    plot((mean(lh(bid))-lhd)./scales(1),(mean(lv(bid))-lvd)./scales(3),'go','markerfacecolor','g');
    plot((mean(lh(cid))-lhd)./scales(1),(mean(lv(cid))-lvd)./scales(3),'go','markerfacecolor','g');
    plot((mean(lh(did))-lhd)./scales(1),(mean(lv(did))-lvd)./scales(3),'go','markerfacecolor','g');
    plot((mean(lh(zid))-lhd)./scales(1),(mean(lv(zid))-lvd)./scales(3),'go','markerfacecolor','g');
    if gain(2,1) < 0.2 || gain(2,3) < 0.2

    text((mean(lh(aid))-lhd)./scales(1),(mean(lv(aid))-lvd)./scales(3),'-1,0');
    text((mean(lh(bid))-lhd)./scales(1),(mean(lv(bid))-lvd)./scales(3),'1,0');
    text((mean(lh(cid))-lhd)./scales(1),(mean(lv(cid))-lvd)./scales(3),'0,-1');
    text((mean(lh(did))-lhd)./scales(1),(mean(lv(did))-lvd)./scales(3),'0,1');
    text((mean(lh(zid))-lhd)./scales(1),(mean(lv(zid))-lvd)./scales(3),'0,0');
    end

    plot((mean(rh(aid))-rhd)./scales(2),(mean(rv(aid))-rvd)./scales(4),'ro','markerfacecolor','r');
    plot((mean(rh(bid))-rhd)./scales(2),(mean(rv(bid))-rvd)./scales(4),'ro','markerfacecolor','r');
    plot((mean(rh(cid))-rhd)./scales(2),(mean(rv(cid))-rvd)./scales(4),'ro','markerfacecolor','r');
    plot((mean(rh(did))-rhd)./scales(2),(mean(rv(did))-rvd)./scales(4),'ro','markerfacecolor','r');
    plot((mean(rh(zid))-rhd)./scales(2),(mean(rv(zid))-rvd)./scales(4),'ro','markerfacecolor','r');
elseif plottype == PLOTSD
    plot(lh(good)-fx(good),stds(good,1),'go');
    hold on;
    plot(rh(good)-fx(good),stds(good,2),'ro');
    plot(lv(good)-fy(good),stds(good,3),'go');
    plot(rv(good)-fy(good),stds(good,4),'ro');
    legend('L','R');
elseif plottype == VERT
    plot(rh,'.');
    hold on;
    if threeeyes
    plot(xh,'r.');
    else
        
    plot(lh,'r.');
    end
else
plot(mean(fx(aid)),median(rh(aid)),'.');
hold on;
plot(mean(fx(bid)),median(rh(bid)),'.');
if threeeyes
plot(mean(fx(aid)),median(xh(aid)),'r.');
plot(mean(fx(bid)),median(xh(bid)),'r.');
end
end

if calcnoise
    noise = abs(median(sds)./scales);
    fprintf('rms noise LH %.4f, RH %.4f,LV %.4f, RV %.4f,XH %.4f, XV %.4f\n',...
        noise(1),noise(2),noise(3),noise(4),noise(5),noise(6));
end

if mkcal
    scales(1,1) = 1./gain(2,1); %lh
    scales(1,2) = 1./gain(1,1); %rh
    scales(1,3) = 1./gain(2,4); %lv
    scales(1,4) = 1./gain(1,4); %rv
    scales(1,5) = 1./gain(3,1);
    scales(1,6) = 1./gain(3,4);
    scales(2,1) = mean(lh(allid)-fx(allid));
    scales(2,3) = mean(lv(allid)-fy(allid));
    scales(2,2)  = mean(rh(allid)-fx(allid));
    scales(2,4) = mean(rv(allid)-fy(allid));
    scales(2,6) = median(xv(allid)-fy(allid));
    scales(2,5) = median(xh(allid)-fx(allid));
    gain = scales;
end

if summarize
    pos.Data = pos;
    pos.Stds = stds;
    pos.ChanLabel = {'Left H'  'Right H'  'Left V'  'Right V'  'Fix X'  'Fix Y'};
end

function HitPoint(a,b,id)

fprintf('%d\n',id);
