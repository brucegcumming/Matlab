function ScrollSpikes(src, evnt)DATA = GetDataFromFig(src);if src ~= gcf    return;endDATA = GetDataFromFig(src);if isempty(DATA.spklst)    DATA.spklst = 1:100;endif DATA.plotspk.bytrial    nt = DATA.currenttrial + evnt.VerticalScrollCount;    DATA.currenttrial = AllV.PlotTrialSpikes(DATA,nt);    set(DATA.toplevel,'UserData',DATA);    return;elsespk = minmax(DATA.spklst);w = DATA.spksperview;if sign(evnt.VerticalScrollCount) > 0    spk(1) = spk(2);    spk(2) = spk(2)+w;elseif sign(evnt.VerticalScrollCount) < 0     spk(2) = spk(1);    spk(1) = spk(1)-w;endendif spk(2) > 0 && spk(1) < DATA.neventsDATA.spklst = spk(1):spk(2);AllV.PlotSpikes(DATA,DATA.spklst);endset(DATA.toplevel,'UserData',DATA);    