classdef PC
properties
        CurrentVersion = '1'
        STIM_GRATING = 3;
end
    methods (Static)        
 varargout = PlotClusters(name, varargin)

 out = doentry(DATA, Clusters, id)
FitDriftMeanSpike(M, varargin)

 C = UpdateWithNewCut(C, DATA)
 DATA = SetDefaults(DATA)

 X = PopupWindow(a, b, fcn, varargin)
 TagMenu(a, b, fcn)
 Expt = combine(DATA, varargin);
 
 GuiMenu(a,b, type)

 DATA = AutoCompare(DATA, varargin)
 
 SetGUI(DATA)

 [C, ok] = GetClusterInfo(C, cl, varargin)

 DATA = CalcDistances(DATA)

 PlotDuplicates(DATA, varargin);
 [d, drop, C] = DistanceMeasure(C, cl, type, varargin)

X = CompareClusters(Ca, Cb, DATA, varargin)

GetComment(DATA,a,b)

 AddComment(a,b,str)

 PrintComments(DATA,e,p)

 HitTrial(data,b, cell)

 PlotISIHist(DATA, c);
 
 SetRateseqPlot(a,b,cellid)

 PlotMenu(a, b, fcn, type, varargin)

 HitShapePlot(a,b,id, c)

 PlotClusterRates(DATA, type,varargin)

 DATA = CalcCellMeans(DATA)

 E = CountSpikes(Expt, C, clnum,xcl)

 PlotCellRates(DATA,type)

 MarkExpts(DATA,type)

 dup = isduplicate(DATA, row, p, cl)

 [true, cellid, clid] = isacell(DATA, row, p, clid)

 SavePlotClustersConfig(DATA,file, varargin);

 DATA = ApplyConfig(DATA, varargin)

 DATA = ApplyLayout(DATA, varargin)

 OptionMenu(a, b, fcn)

out = ReplotXcorrs(a,b, type)

  SetCheckExclusive(a)

 xc = ShapeCorr(P,Q, varargin)

 PlotExpts(DATA)

 FindDuplicates(DATA, cell)

 PlotXcorrs(DATA, xcorrs, expts, bycell)

 str = ProbeLabel(p, DATA, cl)

 args = PlotArgs(DATA)

 Update(a,b)

 [str, value, it] = GetPopString(tag, varargin)

 [value, it] = GetCheck(tag, varargin)

 id = FindSpikes(DATA,C, xcl)

 FixBoundary(DATA, Clusters, e, p)

 plots = PlotClusterPoints(C, uid, cid, varargin)

  h = AddTitle(DATA, C, titlemode)

 FastAxes(ax)

 plots = PlotClusterXY(DATA, C, varargin)

 [r,y] = CalcRadius(E,xy)

 x = Rprime(r)

 xy = XYSpace(C)

 [x, nsp, crit] = PlotHist(xy, varargin)

 Expt = CountExptSpikes(DATA,Expt,C,clnum)

 Expt = PlotExptCounts(DATA, e, p, cl, varargin)

  C = GetCurrentCluster(DATA)

  DATA = ConditionalPlotXY(DATA, C, force, varargin)     

 DATA  = ShowData(DATA, ex,p, varargin)

 result = PlotClusterHistogram(DATA, C, refit, varargin)

 AddTrigHist(DATA, C, cl, varargin)

 PlotCorrelogram(C, varargin)

 h = PlotMeanSpike(C, p, cluster, varargin)

 HitImage(src,b, type)

 Expt = PlotCombinedExpt(DATA)

 [Expt, AllExpt] = PlotSelectedExpts(DATA, varargin)

 CellChanged(DATA)

 PlotXcorr(a,b, pa, pb)

 Q = Cell2Cluster(cell,Clusters)

 HitXcorrAll(a,b, type, cells, expts)

 HitXcorr(a,b, id, ex, cells)

 synci = SyncIndices(xc)

 NewCellSummary(a,b,e,p)

 PlotCellSummary(DATA, e, p)

 [xc, synci, allxc] = meanccf(DATA, id, a,b)

 d = meanmahal(DATA, id, a)

 ms = MeanSpike(DATA, id, a)

 HitPoint(a,b, C, mode, pt)

 HitPopPoint(a,b, ex, p, cell)

 Clusters = ReloadClusters(DATA, eid)

 rawxy = RawXY(C, xy)

 [Clusters, details] = LoadClusters(name)

 DATA = LoadTrialLists(DATA)

 DATA = LoadCellFile(DATA)

 nt = PlotTrialSpikes(DATA, nt, varargin)

 SelectTrial(src, b, type)

 PlotTimeRanges(DATA)

 ChangeTag(a,b,fcn)

 [f, isnew] = SetFigure(DATA, tag, varargin)

 DATA = ReFitAll(DATA, fittype)

 ReFit1D(DATA, ratio)

 DATA = ReFit3means(DATA, ratio)

 DATA = ReFitGMDip(DATA, varargin)

 DATA = LoadExtra(DATA, force)

 [CC, GM] = CondenseClusters(C)

 [CC, GM] = CondenseCluster(C)

 SaveExtras(DATA)

 [dp, fits, details] = Fit2Gauss(C, r, varargin)

 [theta, c, details] = BestAngleGM(xy, G, dipfit, varargin)

 tag = GetFigureTag(src)

 XYCluster(src,b, type, varargin)

 DATA = SaveCluster(DATA, pt, quick, varargin)

 in = InGraph(pt, ax)

 xy = AxPos(ax, pos)

 XYButtonPressed(src, data)

 distance = DistanceToEllipse(E, pos);

 XYButtonReleased(src, data)

 C = ClassifySpikes(DATA, E, mode)

 ScrollWheel(src, evnt)

 FinishXYPlot(ax, DATA, e,p)

 XYButtonDragged(src, data)

 CellTrackMenu(a,b, type)

 DATA = ExcludeTrialsForCell(DATA, probe, cluster, varargin)

 [e,p] = cPoint(DATA)

 SetTrialList(DATA, C, strial)

 PlotAllCellMean(DATA, type, varargin)

 DATA = SetAllXYClustering(DATA, onoff) 

 MarkTrialStarts(Expt, ticks, xcl)

 sdx = PlotXYSequence(DATA, pos, varargin)

 HitXYseq(a,b)

 [eid, cid, clid] = FindCell(DATA, cellid, expts)

 PlotExptCells(DATA, type)

 varargout = PlotAllCell(DATA, type, varargin)

 CallAllVPcs(DATA, eid, pid, varargin)

 SetCellFromLine(a,b, cluster, cell, varargin)

 DeleteCellFromLine(a,b, cluster, cell)

 SetCellFromSubplot(a,b, cell)

 DATA = PlotExptsProbe(DATA, type);

 [S, new] = CheckAllSpikes(DATA, e,p, varargin)

 need = ClusterNeedsRefresh(C, I)

 [Clusters, DATA, e] = CheckClusterLoaded(DATA,e, varargin)

 S = AddSpikes(DATA, S, e,p, varargin)

 PlotAllExpt(DATA, type)

 PlotAllProbeMean(DATA, type, varargin)

 PlotAllExptProbeMean(DATA, type, varargin)

 cmenu = AddCellContextMenu(DATA, type)

 cmenu = AddContextMenu(DATA, type)

 cmenu = AddLineContextMenu(DATA, h, e, p, varargin)

 h = PlotProbeMeans(C,type, varargin)

 DATA = PlotAllProbe(DATA, type)

 mspid = GetMeanSpikeProbe(C, p)

 plotpos = SetPlotPos(DATA,np, nr, nc)

 ReplotXY(a,b,eid, probe,cid)

 DATA = PlotAllProbeXY(DATA,varargin)

 handles = AddCellLabels(DATA, eid, probe, varargin)

 missed = MissedCell(DATA, pos)    

 HitExptPlot(src, b, type, e)

 HitXYPlot(src, b, e,p)

 TightPlot(ax)

 PlotAllCellXY(DATA)

 h= DrawLine(E,varargin)

 [nx, ny] = Data2Norm(x,y)

 h= testDrawLine(E,varargin)

 h= testDrawEllipse(E,varargin)

 h= DrawEllipse(E,varargin)

 PlotAllCellSpikes(DATA)

 DATA = SetCellEntry(DATA, C,  e, p, c, cellid, varargin)    

 SetCellNumber(a,b, fcn)
cells = ClusterList(DATA, Clusters)
 SaveCellList(DATA, varargin)

 DATA = PlotCellList(DATA, varargin)

 DATA = InitInterface(DATA)

 AddPlotMenu(sm, type)

 AddExptList(hm, callbacklabel, DATA)

 DATA = LoadXcorrFiles(DATA)    

 DATA = LoadFullVs(a,b, type, varargin)

 DATA = LoadAll(a,b, type, varargin)

 [Expt, DATA]  = LoadExpt(DATA, e)    

 [DATA, Expts]  = LoadExpts(DATA, varargin)

 Clusters = FixClusters(Clusters)

 DATA = CheckExpts(DATA, type)

 res = CheckSpikeFiles(DATA, type)

 DATA = CheckClusterLineSign(DATA)

 X = CheckClusters(DATA, type)

 DATA = ConvertExclusion(DATA)

 DATA = CheckExclusion(DATA, varargin)

 name = FileName(DATA, ex, probe, type)

 DATA = CompareShapes(DATA, type)    

 SpikeDistances(DATA, type)

 [sizes, d, quality] = SpikeDistance(DATA, eid)

 CompareProbes(DATA, type)

 xcs = ShiftXcorr(allshape, a, b, npts)

 DATA = PlotShapes(DATA, type)

 PlotShape(DATA, spk)


 result = CompareProbesShape(DATA, ex, type)  

 CellFinderTest(DATA, type)

 DATA = FindMissing(DATA)

 DATA =  FillCellList(DATA, mode)

 CellFinder(DATA)

 [CellId, details] = MakeCellId(cellps)

 DATA = TrackCluster(DATA, spk, ex, varargin)

 h = DrawBoxes(DATA, imtype, varargin)

 h = DrawBox(ex, p, imtype, varargin)

 xc = CalcTemplateXcorr(tmpl, C)

 EstimateDrift(DATA, where)

 CompareShape(DATA, Clusters, spk)

 DATA = MarkCurrentCluster(DATA)

 DATA = ClearSelections(DATA, force, setcurrent)

 DATA = SetProbe(mn,b, dir)

 DATA = SetExpt(mn,b,  fcn)

 DATA = NextButton(mn,b, dir)

 DATA = SpoolCurrentSpikes(mn,b)

 res = SpoolSpikeFile(DATA,e,p)

 [x, details] = FindExcludedTrials(DATA,e,p, cluster, C)

 go = CheckSpikeFile(DATA, C)

 Get2DMaxima(DATA)

 name = SpikeFileName(DATA, C)

 [handles, details] = PlotPopPoints(X,Y, varargin)

 PlotMahalImage(DATA, type, varargin)

 RateMenu(src, ks, fcn)

 RateSeqKeyPressed(src, ks, fcn)

 DATA = RateZoom(DATA,inout, cell)

 [yl, details] = GetYRange(h, xl)

 ch = AddClusterString(DATA, h, cl)

 XYKeyPressed(src, ks, fcn)

 DATA = FlipLineCrit(a,b)

 DATA = AddSelectedCells(DATA)

 KeyReleased(src, ks, fcn)

 KeyPressed(src, ks, fcn)

 DATA = AddKey(DATA, tag, varargin)

 DATA = DeleteCell(DATA, e, p, cl)

 s = OldExLabel(DATA, j)

 s = ExLabel(DATA, j)

 DATA = PlotAllClusters(mn,b, varargin)

 DATA = ReadTemplateResults(DATA, nc)

 bad = CheckFitDim(C)

 [DATA, Clusters] = ReadClusterResults(DATA, Clusters)

 DATA = AddError(DATA, varargin)

 Expts = SetExptTimeOffset(Expts)

 DATA = ReadClustersIntoData(DATA, Clusters, exid)

 CloseTag(tag)

 LoadAllSpikes(DATA, varargin)

 LoadSelectedSpikes(DATA,eid,pid)

 [h, details] = QuickSpikes(DATA, pos, varargin)

 stopped = SpoolAllProbes(DATA, e, Spikes, Clusters, varargin)

 PlotSyncSpikes(DATA, eid, probes, clnum, varargin)

 [voffset, ylim] = SetVOffset(DATA, AllSpikes, e, p)

 [stopped, h] = SpoolSpikes(DATA, pos, varargin)

 ploth = PlotSpikes(DATA, pos, spkid, Spks, C, varargin)

 CheckRates(DATA, expname, cell)

 R = CheckAllRateSequences(DATA)

 RatePlotMenu(a,b,fcn)

 PlotRateCheck(R)

 HitExptPoint(a,b, ex, cell)

 res = CheckExptsRateSequence(Expts)

 [err, counts] = CheckExptRateSequence(Expt)
DATA = CompareAutoClusters(DATA, expts)

PlotAllClusterMean(DATA, type, varargin)

    end
end
