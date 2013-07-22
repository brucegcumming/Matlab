function CompareSpikes(varargin)
%If combine and AllVpcs are both up, compares spikes in the two.
showse = 0;
cfig = 1;
cls = 1;
holdon = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'showsem')
        showse =1;
    end
    j = j+1;
end
    
DATA = get(cfig,'UserData');
fullvname = sprintf('X:/Utah/jbe/G040/Expt24.p%dFullV.mat',DATA.probe);
load(fullvname);
cls = unique(DATA.AllData.Spikes.codes);
nevspks = length(DATA.AllData.Spikes.codes);
FullV.t = AddFullVIndex(FullV);
GetFigure('CompareSpikes');
if ~holdon
    hold off;
end
colors = mycolors;
for j = 1:length(cls)
    cl = cls(j);
sid = find(DATA.AllData.Spikes.codes == cl);
t = DATA.AllData.Spikes.nevtimes(sid);
tid = find(ismember(FullV.t,t));
spts = 1:30;
allpts = repmat(spts',1,length(tid))+repmat(tid,length(spts),1);
AllV = FullV.V(allpts);
meanvs{j} = mean(AllV,2); 
svars{j} = var(double(AllV),[],2);
plot(meanvs{j},'color',colors{j});
hold on;
scale = DATA.AllData.Spikes.spkscale.*10;
meanv{j} = mean(DATA.AllData.Spikes.values(sid,:),1); 
vars{j} = var(DATA.AllData.Spikes.values(sid,:),[],1);
plot(meanv{j}.*scale,'--','color',colors{j});
tids{j} = tid;
if showse
    h = errorbar(1:length(meanv{j}),meanv{j}.*scale,sqrt(vars{j}).*scale);
    set(h,'color',colors{j},'linestyle','--');
end
end
dps = (meanvs{1}-meanvs{2})./sqrt(svars{1}+svars{2});
dp = (meanv{1}-meanv{2})./sqrt(vars{1}+vars{2});
tid = tids{2};
title(sprintf('DP%.1f',sum(abs(dp))));
GetFigure('CompareSpikeHist');
if ~holdon
    hold off;
end
spk = meanvs{2}-mean(meanvs{2});
x = conv(spk,double(FullV.V));
y = diff(sign(diff(x)));
nev = length(tid) * 3;
maxid = find(y < 0);  %maxima for template score
[a,b] = sort(x(maxid),'descend');
[y,x] = hist(a(1:nev),500);
hold off;
bar(x,y);
y = histc(double(AllV)' * spk,x);
hold on;
plot(x,y,'r');
title(sprintf('%d/%d Nev. Showing %d Events\n',length(tid),nevspks,nev));
%[y,x] = hist(double(AllV)' * spk,500);
%plot(x,y,'g');
GetFigure('CompareSpikeTimes');
if ~holdon
    hold off;
end
e = DATA.exabsid(1);
p = DATA.probe;
cx = DATA.AllClusters{e}(p).cx(sid);
cx = (cx-mean(cx))./std(cx);
plot(DATA.AllData.Spikes.times(sid,:),cx,'b.','buttondownfcn',@HitPoint);
hold on;
cname = [DATA.datadir '/Expt' num2str(DATA.exabsid(1)) 'ClusterTimes.mat'];
dname = [DATA.datadir '/Expt' num2str(DATA.exabsid(1)) 'ClusterTimesDetails.mat'];
if exist(cname)
    load(cname);
    if exist(dname)
        load(dname);
        Ts = DATA.Expts{e}.Trials;
        for j = 1:length(Ts)
            plot([Ts(j).Start(1) Ts(j).End(end)],[0 0],'k-');
        end
        id = find(ClusterDetails{p}.clst == cl+1);
        xy = ClusterDetails{p}.xy(id,1);
        x = (xy(:,1)-mean(xy(:,1)))./std(xy(:,1));
        plot(ClusterDetails{p}.t(id).*10000,x,'r.','buttondownfcn',@HitPoint);
        nspk = length(x);
    else
        plot(Clusters{DATA.probe}.times.*10000,1.9,'r.');
        nspk = length(Clusters{DATA.probe}.times);
    end
end
title(sprintf('%d Nev Spikes, %d AllVpcs',length(sid),nspk));


function HitPoint(a,b)

a = get(gca,'currentpoint');
t = a(1,1);
fprintf('Spike at %.4f\n',t./10000);
AllVPcs('plottrialtime',t);

function t = AddFullVIndex(FullV)
first = 1;
vt(sum(length(FullV.blklen))) = 0;
for j = 1:length(FullV.blklen)
    offset = round((FullV.blkstart(j)-FullV.blkstart(1)).*30000.237);
    last = first+FullV.blklen(j)-1;
    vt(first:last) = offset+[1:FullV.blklen(j)];
    first = last+1;
end
t = vt;

function t = AddFullVTimes(FullV)
if isfield(FullV,'samper')
    FullV.interval = FullV.samper;
end
FullV.interval = 1./30000.237;
    first = 1;
    vt(sum(length(FullV.blklen))) = 0;
    for j = 1:length(FullV.blklen)
        last = first+FullV.blklen(j)-1;
        vt(first:last) = FullV.blkstart(j)+[1:FullV.blklen(j)].*FullV.interval;
        first = last+1;
    end
    if length(vt) > size(FullV.V,2)
        PrintMsg(DATA.logfid,'sum blklen > length(V)\n');
        vt = vt(1:size(Vall.V,2));
    end
    t = vt;