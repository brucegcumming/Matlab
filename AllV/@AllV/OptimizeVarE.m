function C = OptimizeVarE(DATA)       AllV.SetFigure(DATA.tag.vare, DATA);   subplot(2,2,1);    x = DATA.energy(1,:);    y = DATA.spkvar(DATA.probe(1),:)./DATA.energy(1,:);    Cs(1) = AllV.CutAndPlot(x,y,DATA.energy(1,:));    drawnow;    subplot(2,2,2);    x = DATA.energy(1,:).^2;    y = DATA.spkvar(DATA.probe(1),:).^2 ./DATA.energy(1,:).^2;    Cs(2) = AllV.CutAndPlot(x,y,DATA.energy(1,:));    drawnow;    C = Cs(1);    C.sign = 0;