function aid = PlotArtifacts(DATA)

aid = [];    
Spks = DATA.AllData.Spikes;
id = find(Spks.codes(:,2) == 8)
amean = mean(Spks.values(id,:));
scores = Spks.values * (amean-mean(amean))';
GetFigure('Artifacts');
hold off; 
x = max(Spks.values(:,10:15),[],2);
plot(x,scores,'.','markersize',1);
a = questdlg('Classify?','popup','Cancel','OK','OK');
if strcmp(a,'OK')
E = AddEllipse(gcf,'wait','color','r','line','timeout',30);
if ~isempty(E)
C = XYClassify(x,scores,E);
hold off; 
plot(Spks.values(C.id,:)','color','b');
hold on;
plot(amean,'r','linewidth',2);
plot(mean(Spks.values(C.id,:)),'g','linewidth',2);

a = questdlg('Apply?','popup','Cancel','OK','OK');
if strcmp(a,'OK')
aid = C.id;
end
end
end


