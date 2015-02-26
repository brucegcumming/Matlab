function t = CompareNevNsx(nev, nsx, ids)
cfig = 1;
probe = 4;
coff = 2;

colors = mycolors;
pid = find(nev.Data.Spikes.Electrode == probe);
if nargin > 2
    pid = pid(ids);
end
t = double(nev.Data.Spikes.TimeStamp(pid));
spts = -30:50;
allpts = repmat(spts',1,length(t))+repmat(t',length(spts),1);
AllV = nsx.Data(allpts);
GetFigure('NevNsx');
plot(spts,mean(AllV,2),'color',colors{1+coff});
hold on;
plot(mean(nev.Data.Spikes.Waveform(pid,:)),'color',colors{2+coff});
%plot(AllV);


function t = AddFullVTimes(FullV)
if isfield(FullV,'samper')
    FullV.interval = FullV.samper;
end
FullV.interval = 1./30000;
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