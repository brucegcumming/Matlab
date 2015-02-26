function DATA = SpkVMenu(a,b, type)
DATA = GetDataFromFig(a);

if type == 1
a = max(DATA.AllData.Spikes.values(DATA.spklist,:),[],2);
vm = min([prctile(a,99.5) .* 1.1 5]);
DATA.plot.SpikeMaxV = vm;
a = min(DATA.AllData.Spikes.values(DATA.spklist,:),[],2);
vm = max([prctile(a,0.5) .* 1.1 -5]);
DATA.plot.SpikeMinV = vm;
if isempty(b) %called from gui, so fix the ranges
DATA.plot.autoVrange = 0;
end
set(DATA.toplevel,'UserData',DATA);
elseif type ==2
Spks = GetSpikeStruct(DATA);
id = find(Spks.codes(:,2) < 8);
a = cmb.FitGaussMeans([DATA.Spikes.cx(DATA.spklist)' DATA.Spikes.cy(DATA.spklist)'],2, 'verbose');
DATA.Expts{DATA.currentexpt(1)}.Clusters{1,DATA.probe}.mahal = a.mahal;
C = DATA.Expts{DATA.currentexpt(1)}.Clusters{1,DATA.probe};

hold on;
ezcontour(@(x,y)pdf(a.obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));
plot(a.obj.mu(1,1),a.obj.mu(1,2),'+');
plot(a.obj.mu(2,1),a.obj.mu(2,2),'+');

if isfield(Spks,'pcs')
a = cmb.FitGaussMeans(Spks.pcs(DATA.spklist,:),2);
fprintf('Distance for %d PCs %.2f\n',size(Spks.pcs,2),a.mahal);
C.pcmahal = a.mahal;
[a,b,c] = cmb.BestAngle(Spks.pcs(DATA.spklist,1), Spks.pcs(DATA.spklist,2),3);
fprintf('PC1/2 Bestangle %.1f, %.3f H%.4f\n',a * 180/pi,b,c.hdip);
C.pcdip = c.hdip;
C.pcbmc = b;
end
[a,b,c] = cmb.BestAngle(DATA.Spikes.cx(DATA.spklist)', DATA.Spikes.cy(DATA.spklist)',3);
fprintf('Bestangle %.1f, %.3f H%.4f\n',a * 180/pi,b,c.hdip);
C.dip = c.hdip;
C.bmc = b;
DATA.Expts{DATA.currentexpt(1)}.Clusters{1,DATA.probe} = C;
set(DATA.toplevel,'UserData',DATA);

elseif type ==3  %recalc PCA
DATA = CheckForPCA(DATA,DATA.spklist, 1);
set(DATA.toplevel,'UserData',DATA);
end

