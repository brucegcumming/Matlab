classdef cmb
properties
 Version = '1.1'
end
methods (Static)
varargout = combine(name, varargin)
DATA = Init(DATA, name, TOPTAG, layout)
[out, varargout] = SpkCache(DATA, e,p, varargin);
CloseCombine(a,b, varargin)
layout = ReadLayout(file)   
DATA = SetCluster(DATA, eid, probes, varargin);
DATA = BuildAllTimes(DATA, varargin);
DATA = CheckState(DATA, varargin);
SaveWithEM(a,b)
DATA = CalcOnlineShapes(DATA)
PlotAdjacentXcorrs(DATA, probes, type)
DATA = PlotXcorr(DATA, probes, type)
DATA = ReadGridFile(DATA)
DATA = ReadFile(DATA, name, varargin)
DATA = CheckCellExptList(DATA)
DATA = LoadCellFile(DATA)
probesource = FindProbeSources(DATA)
PrintLogData(DATA, type)
s = ExptComments(Expt)
DATA = PrintComments(DATA)
DATA = AutoCut(DATA, eid, idx, varargin)
NotBusy(DATA)
ShowLFPPlot(DATA)
[figpos, F] = SetFigure(tag, DATA, varargin)
SetCellChooser(DATA)
AllCellKeyPressed(src,ks, op)
P =  NextPlot(P, AllCellRes, step)
AllCellPlots(a,b, op)
str = ShowString(DATA, Expt)
Expts = CountTxtSpikes(Expts, probe, cl)
Spks= GetSpikeStruct(DATA)
DATA = LoadClusters(DATA, cfile, varargin)
DATA = ClearExptClusters(DATA)
DATA = LoadOnlineClusters(DATA, cfile)
[DATA, Stimvals, needv, retval] = CheckCombine(DATA, interactive)
AutoClusters = SaveAutoClusters(DATA)
OptimizeDprimeHit(a,b)
C = OptimizeDprime(DATA)
dprime = ClusterDprime(x, DATA, nspks)
Cluster = SmallCluster(AllClusters, e, p)
[Expt, DATA, plotres] = CombinePlot(DATA, dlgs, varargin)
[Expt, res] = PlotCombined(DATA, Expt, varargin)
SetProbeList(DATA)
DATA = ReadFinishedDir(DATA, name, varargin)
DATA = ReadDir(DATA, name, varargin)  %% Online Plots
PlotSpike(DATA, ispk, probe)
PlayNextSpike(a,b)
PlotISIPair(DATA, pair)
CutTrial(a,b)
PlayLastTrial(a, b)
PlayNextTrial(a, b)
SelectTrial(a, b)
PlayOneTrial(DATA, a, b, varargin)
PlotLFPRaw(state, Trial, crrate);
DATA = PlotTrialSpikes(DATA, itrial, colors, clusters)
DATA = PlotSpikes(DATA, ispk)
DrawSpikeWaves(DATA, ispk, nclusters, ctype)
[ispk, sspk,cx, cy] = PlotTrialSyncSpikes(DATA, times, probes, colors, varargin)
[DATA, ispk] = APlotTrialSpikes(DATA, times, colors, nclusters, classify, varargin)
xc = ExtraLabels(Trial)
spikelist = WhichClusters(top,varargin)
[DATA, counts] = CountSpikesB(DATA, expid, pid, spikelist, varargin)
[DATA, counts] = CountSpikes(DATA, expid, varargin)
[eid, cid, shifts] = FindNoCell(DATA, probe, expts)
[eid, cid, clid] = FindCell(DATA, cellid)
nc = MaxSpikeCode(DATA, expspks)
nc = CountClusters(Cluster)
outname = Expt2FileName(DATA,Expt, cluster)
[outname, prefix] = CombinedName(DATA,eid,cluster, varargin)
savecl = SetExptClusters(caller,b, varargin)
DelClusterButton(caller,b)
ClrSpkWin(caller,b)
cfile = ClusterFile(DATA,varargin)
cfile = CombinerLst(DATA)
[x,y, DATA] = GetClusterSpace(DATA, Expt)
DATA = DrawClusters(DATA, cluster, setfig)
CheckSpoolButton(DATA)
cluster = CopyClusters(cluster, ec)
h =  DrawCluster(cluster, color)
DATA = MainMenu(a,b,type)
DATA = OptionMenu(a,b,type)
DATA = SpkVMenu(a,b, type)
[theta, c, details] = BestAngle(x,y, test, varargin)
c = BimodalCoeff(x, e)
res = FitGaussMeans(X,N, varargin)
SetTrialRange(a,b, type)
SetXYCluster(a,b, type, plot, val,varargin)
PlotEig(DATA)
DATA = Plot3DClusters(a,b);
AddClusterIdButton(a,b)
DATA = PlaySpikes(DATA, expid, varargin)
PlotPCAs(DATA, eid)
OUT = CalcTrodeCorrs(DATA,eid)
PlotTrodeXcorr(a,b)
SetSetBox(DATA, ei)
RescaleClusterPlot(a,b)
SpoolSpikes(a,b)
NextList(a,b, varargin)
ClearMouse()
PlotXYDensity(energy,vw)
DensityPlot(a,b)
DATA = cListExpts(DATA, Expts, varargin)
[id, ok] = CheckListForName(list,name)
CheckLists(DATA)
cplists(caller,b)
DATA = ListSubExpts(DATA, id, varargin)
timerfna(tim, varargin)
timerfn(tim, varargin)
LoadAllProbes(a,b, varargin)
DATA = GetAllProbeFig(DATA)
PlotAllExptsAndProbes(a,b, varargin)
PlotSubspace(a,b, varargin)
SpoolAllExpts(a,b, toplevel, type)
HitXcorrPt(src, b, id, j)
HitXcorr(src, b, id)
PredictChoices(DATA);
PlotXcorrs(DATA, plottype)
PlotRateSequences(DATA, varargin)
PlotAllCells(a,b, plottype)
PlotAllExpts(a,b, varargin)
MarkProbes(a,b)
HideSpikes(a,b)
LoadAllProbeParams(a,b, varargin)
DATA = LoadAllProbesSpikes(DATA)
DATA = LoadAllProbeSpikes(a,b, varargin)
LoadSpikeTimes(a,b, varargin)
xc = CalcXcorrDC(DATA, eids, sa, sb)
DATA = ReClassifyAll(DATA, varargin)
xc = CalcXcorrV(DATA, eids, sa, sb)
xc = CalcXcorr(DATA, eids, sa, sb)
TrackTemplates(DATA)
lProbeClusters(a,b)
DATA = PlotAllProbeXY(a,b)
DATA = SetXYRanges(DATA, cx, cy)
[xr, yr] = ClusterRange(C, p)
DATA = AddMultiProbeGUI(DATA)
ListCells(a,b)
MakeCellList(a,b, varargin)
SetCellNum(a,b)
MakeCellTemplate(a,b)
RePlotCellList(a,b)
AddOneCellToList(caller,b,varargin)
AutoFillCellList(caller,b, varargin)
[cells, probes] = CellsSelected(DATA)
AddCellToList(a,b, varargin)
[a,b] = TrialRange(DATA)
SaveCellList(a,b, varargin)
ClearAnnotations(a)
PlotNewCellList(DATA, varargin)
DATA = MarkCurrentCell(DATA)
HitImage(a,b, type)
str = ProbeLabel(Expt)
PlotCellList(DATA, varargin)
AddProbeList(DATA)
GuiResize(a,b);
DATA = BuildGUI(DATA)
DATA = DoConfig(a,b, fcn)
DATA = DoLayout(a,b, fcn)
MakeTemplates(a,b, varargin)
PlotTemplates(a,b, varargin)
AddTemplate(a,b, varargin)
GetCellNumber(DATA, eid, probe)
ReloadClusters(a,b, varargin)
SetCellByTemplateScore(DATA, cell, probes, trials, tn)
CalcSpikeShapes(a,b, varargin)
PlotTemplateScores(DATA,ti,varargin)
DATA = SaveSpikeShape(DATA, outname)
AddComment(a,b)
DATA = CalcMeanSpike(DATA,expid)
Spks = GetSpikeData(DATA, p)
ReCombineByName(DATA, name, backup)
out = ReCombineAll(a,b, varargin)
spkstats = GetSpkStats(DATA)
AllExpt = CombineAllCells(a,DATA, varargin)
cp = ExptCP(Expt)
CombineAll(a,DATA, varargin)
FitButton(a,b)
DATA = AddFitToData(DATA, plotres, fit)
cntrl_box = setshow(DATA, tag)
str = vec2str(x)
SpikeUIMenu(a,b);
cntrl_box = setoptions(DATA, tag)
ShowUpdate(a,b, varargin)
[value, it] = GetShow(DATA,tag)
SetClusterCheck(DATA)
SetGui(DATA);
SetClusters(a,b,tag)
xcorrhit(a,b, type, varargin);
SetProbeHit(a,b, varargin)
[DATA, ok] = SetProbe(DATA, probe)
[DATA, D] = CheckClusterLoaded(DATA, eid, pid, varargin)
[DATA, D] = ReadCluster(DATA, eid, pid, varargin)
[dips, diperrs] = GetDips(C)
filename = GetProbeFilename(DATA, eid, probe)
DATA = LoadSpikes(DATA, eid, varargin)
DATA = GetNS5Spikes(DATA, eid, probe);
DATA = SetSpkLists(DATA)
Spikes = GetProbeFiles(DATA, probe, subprobe, varargin)
Spikes = GetProbeSpikes(All, filename, varname, probe)
Update(a,b)
UpdateItem(a,b,varargin);
SetField(parent, tag, value, varargin)
value = GetField(tag, varargin)
[value, it] = GetCheck(tag, varargin)
[success] = SetCheck(tag, value,  varargin)
args = PlotArgs(DATA, Expt,varargin)
mousept = myellipse(mousept, finish)
ScrollWheel(src, evnt)
ScrollTrial(src, evnt)
KeyPressed(a, ks)
KeyReleased(a, ks)
NewCluster(a)
DeleteCluster(cl,callfig, varargin)
FigButtonReleased(src, data)
ExcludeTrials(DATA, mousept)
mousept = myrect(mousept, finish)
FigButtonPressed(src, data)
FigButtonDragged(src, data)
ButtonPressed(src, data)
distance = DistanceToCluster(C, pos);
PlotSpikeXY(DATA, spkid, color)
PlotSpikeXYClusters(DATA, expspks, clusters, varargin)
ShowCellLabels(DATA)
PlotDprime(DATA)
ClassifySpikes(mousept,src,varargin)
plotISI(DATA)
[dprime, details] = CalcDprime(DATA,cluster, varargin)
dp = CalcClusterdp(DATA, cl)  %old method using ROC of radial distances.
DATA = DrawXYPlot(DATA, expspks, varargin)
PlotDDF(DATA)
DATA = FinishXYPlot(DATA)
res = isExptCluster(E,c,p)
res = iscluster(C,c,p)
missedtrials = FindMissingTrials(DATA, Vall, ei)
ptsize = CheckPtSize(DATA, nspk)
ButtonReleased(src, data)
ButtonDragged(src, data)
Expts = CheckSaccades(Expts, emfile)
name = BuildName(name)
[aid, bid, n] = FindSync(at,bt,dt,varargin)
[E,V, pc] = SpikePCA(spk, probes, ids)
aid = PlotArtifacts(DATA)
DATA = UpdateAllClusters(DATA)

end
end