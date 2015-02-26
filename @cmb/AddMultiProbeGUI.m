function DATA = AddMultiProbeGUI(DATA)

if ~isfield(DATA,'toplevel') %no GUI 
return;
end
ncw=1./45;
cmb.SetProbeList(DATA);

if ~isempty(findobj(DATA.toplevel,'Tag','MultiProbeMenu'))
return;
end
DATA.show.ed = 1;
if 0
cw = DATA.plot.cw;
ch = DATA.plot.ch;

SPACE = cw;

bp(1) = SPACE;
bp(2) = DATA.gui.toprow+ch+SPACE;
bp(3) = cw.*10;
bp(4) = ch * 1.5;
uicontrol(DATA.toplevel,'Style', 'pushbutton', 'Callback', @cmb.CombineAll, ...
'String', 'Combine All', 'Tag','CombineAll','Position', bp);
bp(1) = bp(1)+bp(3)+SPACE;
uicontrol(DATA.toplevel,'Style', 'pushbutton', 'Callback', @cmb.ReCombineAll, ...
'String', 'ReCombine All', 'Tag','ReCombineAll','Position', bp);
end

it = findobj(DATA.toplevel,'Tag','ProbeId');
bp = get(it,'Position');

if isfield(DATA.probes,'traces')
np = max([DATA.probes.traces]);
DATA.subprobes = np;
if np > 1
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'pop','String',num2str([0:np]'),'Units', 'norm','Position', bp,...
'Tag','SubprobeId','Callback',@cmb.SetProbeHit,'value',1);
%Add extra measures needed for subprobes
[DATA.spkvarnames, DATA.spkvarorder] = GetSpikeVals(DATA,NaN,NaN,NaN,[]);
end
end
bp(1) = bp(1)+bp(3)+0.01;
bp(3) = ncw*3;
uicontrol(gcf,'style','pushbutton','string','>>','units','norm','Position',bp, 'Callback',{@cmb.SetProbeHit,'next'});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'style','pushbutton','string','<<','units','norm','Position',bp, 'Callback',{@cmb.SetProbeHit,'prev'});

hm = uimenu(DATA.toplevel,'Label','&Nchannel','Tag','MultiProbeMenu');
if isfield(DATA,'AllClusters')
uimenu(hm,'Label','&Update Clusters','Callback',{@cmb.LoadAllProbeSpikes, 'allselect'});
else
uimenu(hm,'Label','&Load All Spks','Callback',{@cmb.LoadAllProbeSpikes, 'allselect'});
end
uimenu(hm,'Label','&Load Selected Spks','Callback',{@cmb.LoadAllProbeSpikes,'select'});
uimenu(hm,'Label','&Load All Cluster Params','Callback',@cmb.LoadAllProbeParams);
uimenu(hm,'Label','&Load Spike Times only','Callback',@cmb.LoadSpikeTimes);
uimenu(hm,'Label','&Plot All Spks','Callback',@cmb.LoadAllProbes);
uimenu(hm,'Label','&RePlot All Spks','Callback',@cmb.PlotAllProbeXY);
uimenu(hm,'Label','&Combine All','Callback',@cmb.CombineAll);
uimenu(hm,'Label','&Combine All Cells','Callback',{@cmb.CombineAll,'cells'});
uimenu(hm,'Label','&All Cells + MU','Callback',{@cmb.CombineAll,'mucells','save'});
uimenu(hm,'Label','All Cells + MU quick','Callback',{@cmb.CombineAll,'mucells' 'saveonlyall' 'quick'});
uimenu(hm,'Label','Combine All &MU','Callback',{@cmb.CombineAllCells,'muonly'});
uimenu(hm,'Label','Make Combined File for each cell','Callback',{@cmb.CombineAll,'cells','save'});
uimenu(hm,'Label','&Reombine All Cells','Callback',{@cmb.ReCombineAll,'cells'});
uimenu(hm,'Label','&Reombine All Cells + MU','Callback',{@cmb.ReCombineAll,'mucells'});
uimenu(hm,'Label','&Combine All+LFP','Callback',{@cmb.CombineAll, 'lfp'});
uimenu(hm,'Label','&Recombine All','Callback',@cmb.ReCombineAll);
uimenu(hm,'Label','&Recombine All+','Callback',{@cmb.ReCombineAll 'mkall'});
uimenu(hm,'Label','&Recombine All LFPs','Callback',{@cmb.ReCombineAll 'lfponly'});
uimenu(hm,'Label','&Recombine All spks','Callback',{@cmb.ReCombineAll 'spkonly'});
uimenu(hm,'Label','&Recombine One probe','Callback',{@cmb.ReCombineAll 'oneprobe'});
uimenu(hm,'Label','&Recombine Check List','Callback',{@cmb.ReCombineAll 'listonly'});
uimenu(hm,'Label','&Spike Shape','Callback', @cmb.CalcSpikeShapes);
uimenu(hm,'Label','MeanSpike this probe','Callback', {@cmb.combine 'meanspike'});
uimenu(hm,'Label','->One Probe','Callback', {@cmb.SetProbeHit 'ReloadProbe'});
uimenu(hm,'Label','->CellList','Callback', @cmb.MakeCellList);
uimenu(hm,'Label','List by cell','Callback', @cmb.ListCells);
uimenu(hm,'Label','CrossCorr','Callback', {@cmb.xcorrhit, 1});
if isfield(DATA,'AllClusters')
uimenu(hm,'Label','CrossCorr Neighbors','Callback', {@cmb.xcorrhit, 2});
end
hm = uimenu(DATA.toplevel,'Label','&Plot','Tag','PlotMenu');
uimenu(hm,'Label','All Expts - spkxy','Callback',@cmb.PlotAllExpts);
uimenu(hm,'Label','All Expts and probes','Callback',@cmb.PlotAllExptsAndProbes);
uimenu(hm,'Label','Plot Subspace','Callback',@cmb.PlotSubspace);
uimenu(hm,'Label','LFP Plot','Callback',{@cmb.combine, 'maklfp'});
uimenu(hm,'Label','Export Result to workspace','Callback',{@ExportPlots, 'combine'});

sm = uimenu(hm,'Label','&All Expts ');
uimenu(sm,'Label','Rates','Callback',{@cmb.PlotAllCells, 'rates'});
uimenu(sm,'Label','Normalized Rates (max)','Callback',{@cmb.PlotAllCells, 'normmax'});
uimenu(sm,'Label','Normalized &Rates (mean)','Callback',{@cmb.PlotAllCells, 'normmean'});
uimenu(sm,'Label','Normalized image (max)','Callback',{@cmb.PlotAllCells, 'imagemax'});
uimenu(sm,'Label','Normalized &image(mean)','Callback',{@cmb.PlotAllCells, 'imagemean'});
uimenu(sm,'Label','Var','Callback',{@cmb.PlotAllCells, 'var'});
uimenu(sm,'Label','Var image','Callback',{@cmb.PlotAllCells, 'varim'});
uimenu(sm,'Label','Blank','Callback',{@cmb.PlotAllCells, 'blank'});
uimenu(sm,'Label','Blank image','Callback',{@cmb.PlotAllCells, 'blankim'});
uimenu(sm,'Label','Rate Sequence','Callback',{@cmb.PlotAllCells, 'rateseq'});
uimenu(sm,'Label','Rate Sequence (normalized)','Callback',{@cmb.PlotAllCells, 'normrateseq'});
uimenu(sm,'Label','Rate Sequence+offset','Callback',{@cmb.PlotAllCells, 'offsetrateseq'});
uimenu(sm,'Label','SacTrig','Callback',{@cmb.PlotAllCells, 'sacsdf'});
uimenu(sm,'Label','SacTrigIm','Callback',{@cmb.PlotAllCells, 'sacsdfim'});
uimenu(sm,'Label','test','Callback',{@cmb.PlotAllCells, 'test'});
uimenu(sm,'Label','xcorr-distance','Callback',{@cmb.PlotAllCells, 'xcorr-distance'});
uimenu(sm,'Label','xcorr-cp','Callback',{@cmb.PlotAllCells, 'xcorr-cp'});
uimenu(sm,'Label','xcorr-test','Callback',{@cmb.PlotAllCells, 'xcorr-test'});
uimenu(sm,'Label','Load','Callback',{@cmb.PlotAllCells, 'load'});


