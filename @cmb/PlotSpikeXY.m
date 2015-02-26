function PlotSpikeXY(DATA, spkid, color)
ptsize = cmb.CheckPtSize(DATA, length(DATA.spklist)); %don't set diff sizes for diff clusters!
plot(DATA.Spikes.cx(spkid),DATA.Spikes.cy(spkid),'.',...
'color',color,'markersize',ptsize);

