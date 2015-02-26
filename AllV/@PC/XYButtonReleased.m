function XYButtonReleased(src, data)DATA = GetDataFromFig(src);if DATA.elmousept.down == 0     return;endmode = DATA.elmousept.down;DATA.elmousept.mode = mode;start = get(gca,'CurrentPoint');DATA.elmousept.done = 1;p = DATA.elmousept.pos;DATA.elmousept.down = 0;axdata = get(gca,'UserData');if isfield(axdata,'probe')    pi = axdata.probe;    ei = axdata.eid;    DATA.currentpoint = [ei pi];    fprintf('Clustering E%dP%d\n',ei,pi);else    pi = DATA.currentpoint(2);    ei = DATA.currentpoint(1);endif mode == 1DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; elseif mode == 2  %button was pressed inside ellipse, just move itDATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; endxyr = DATA.elmousept.xyr;set(DATA.toplevel,'UserData',DATA);%touch inside ellispe to make cut. If drawing a line, make cut when%releasedif mode == 2  || (DATA.elmousept.shape  ==1 && mode == 1)%PCCluster(DATA,DATA.elmousept,1);       C =  PC.ClassifySpikes(DATA,DATA.elmousept,1);        if (DATA.plot.hist | DATA.plot.refitgm) && DATA.elmousept.shape == 1            [a,b]  = GMDip(C.xy(:,1),0);            [c,d,e] = GMfit(C.xy,2,1,'idlist',C.clst);            title(sprintf('E%dP%d mahal %.2f,%.2f(1/2)\n',ei,pi,b.mahal(b.best),d));            if DATA.plot.hist                PC.SetFigure(DATA,DATA.tag.hist,'front');                PC.PlotHist(C.xy, C);                title(sprintf('E%dP%d mahal %.2f\n',ei,pi,b.mahal(b.best)));            end        elseif DATA.plot.refitgm && DATA.elmousept.shape == 0            [a,b,c] = GMfit(C.xy,2,1,'idlist',C.clst);            tic;            [d,e, details] = PC.BestAngleGM(C.xy,a,[],'quick');            toc            [aa,bb]  = GMDip(details.xy(:,1),0);            title(sprintf('E%dP%d mahal %.2f (%.2f 1d)\n',ei,pi,b,max(bb.mahal)));        end        PC.PlotExptCounts(DATA);endfigure(src);