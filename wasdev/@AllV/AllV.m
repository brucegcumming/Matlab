classdef AllV
properties
        Version = '1'
    end
    methods (Static)        

        [res, varargout] = AllVPcs(V, varargin)
result = PlotTriggers(Vall, varargin);

 probes = GetProbeList(DATA);


 FullVButtonPressed(a,b)
 
 FullVButtonReleased(a,b)
 ScrollFullV(a,b)
 DATA = mysetappdata(DATA, str, val)

 val = mygetappdata(DATA, str)

 [AllVoltages, DATA] = BuildAllV(DATA, id, spts)

 val = myisappdata(DATA,str)

 DATA = CloseLog(DATA)

 C = CheckClusterFields(C, varargin)

 [V, DATA] = ReadSpikeFiles(DATA, name)

 str = IDStr(DATA)

 C = CheckScoreScaling(DATA, C)

 ClusterFromPoints(DATA)

 type = WhichPlotType(C,clnum)

 [DataClusters, FullVData] = LoadDataClusters(DATA)

 [DATA, Vall, ispk, newdata] = SetupVall(DATA,Vall, ispk, newdata)

  [true, each] =  NeedTemplateForCluster(C, all)

 DATA = GetGuiState(DATA, F)

 DATA = ResetDataForNewProbe(DATA)

 DATA = ReadFromLog(DATA)

 sz = memsize(X)

 S = SmallCluster(C)

 [true, cluster] = ClusterIsSet(C, cl)

 [DATA, ClusterDetails] = LoadClusterDetails(DATA, varargin)

 [DATA, DataClusters, success] = LoadTrigTimes(DATA, checktimes, varargin)

 [id,  th, details] = TriggerV(DATA, rV)

 [res, C,D] = AutoCutOne(DATA, exname, ispk, args)

 res = AutoCutAll(ispk, toplevel, Vall,DATA, args)

 first = PrevCluster(C)

 res = QuantifyQuickClusters(DATA, ispk, varargin) 

 need = NeedMore(C, Evec)

 C = CheckForMean(DATA,C)

 DATA = ReClassify(DATA, varargin)

 DATA = CalcTemplateScores(DATA,  varargin)

 DATA = CalcTemplatesFromMean(DATA, MeanSpike, varargin);

 tt = TimeMark(tt, str, show)

 good = GoodCluster(C)

 C = ClusterFromBoundary(E, C)

 E = BoundaryFromCluster(E, C, n)

 DATA = NextPCs(DATA)

 DATA = SetPCs(DATA, replot, reapply)

 DATA = CalcICA(DATA,ns)    

 [C, Evec, pcs, dip, chspk, errs, details] = CalcPCs(DATA, AllVoltages, nprobepc,  varargin)

 OldPlotFullV(V, ispk, DATA)

 distance  = gmdistance(G)

 [theta, c, details] = BestGMAngle(x,y, test, varargin)

 [theta, c, details] = BestAngleGM(xy, G, dipfit, varargin)

 [theta, c, details] = BestAngle(x,y, test, varargin)

 DATA = CheckClusterMarks(Clusters, DATA)

 CheckClusters(Clusters, str, varargin)

 SaveMeanSpikeOnly(DATA, outname) 

 CheckClusterValues(DATA, C)

 [DATA, id] = SaveClusters(DATA, outname,varargin)

 C = StripClusters(Clusters)

 name = SpkFileName(DATA, varargin)

 probe = ProbeNumber(DATA)

 SaveSpikes(DATA, id, name)

 tcut = IsTemplateCut(E)    

 [E, cluster] = CutAndSave(DATA, varargin)

 [distance, obj, xy, details] = TemplateSpace(DATA, varargin)

 DATA = TemplateGMFits(DATA)        

 [distance, obj, xy, details] = BestSpace(DATA, varargin)

 [G, D, all] = GMfit(X, nd, nr, varargin)

 [G, xy] = IterateTemplateFit(DATA, G)

 [xy, cid, details] = ProjectND(DATA, best, obj)

 Spikes = MakeJamesSpikes(DATA)

 [DATA, details]  = JamesAutoCut(DATA, varargin)

 [E, Scores, tbs, xy, details]  = AutoCut(DATA, varargin)

 DATA = AddErr(DATA,varargin)

 C = CutAndPlot(x,y, energy)

 C = OptimizeVarE(DATA)

 [DATA, E] = OptimizeBoundary(DATA);

 [C, details] = OptimizeClusterBoundary(DATA)

 FullV = SetFullVNames(DATA, FullV)

 PlotXY(xy, clst)    

 c = BimodalCoeff(x, e)

 ChangeCell(a,b,p)

 ChangeProbe(a,b,p)

 OptionMenu(a,b, fcn)

 pos = PlaceUi(a, b ,str)

 args = RetriggerDialog(a,b, fcn)

 FitWindow(F)

 SetMenuChecks(hm, S)

 SelectProbe(a,b,p)

 ProbeSelector(DATA)

 ProbeMenu(a,b, fcn)

 res = PlotClusters(a,b,fcn)

 SetCellFromLine(a,b, cluster, cell)

 DATA = PlotCellList(DATA, varargin)

 DATA = SetCellEntry(DATA, C,  e, p, c, cellid, varargin)    

 SaveCellList(DATA)

 AddSelectorContextMenu(DATA, ax, probe)

 DATA = AddAxisContextMenu(DATA, ax)

 cmenu = AddLineContextMenu(DATA, h)

 M = CalcDistanceMatrices(DATA, nc, varargin)

 cid = AssignCluster(DATA, G)

 SetPlot(a,b, fcn)

 GMMmenu(a,b, fcn)

 GMMButton(DATA, E, fcn)

 DATA = RunGMMFit(DATA,C, varargin)

 PCCluster(a,b, fcn)

 C = GetSubCluster(Cluster, c)

 DATA = CheckTemplates(DATA, C)

 [Expt, matfile] = LoadExpt(DATA, ei)

 Expt = LoadExptA(DATA, exfile, ei)

 res = FitGaussMeans(X,N, varargin)

 [d, details]  = gmdprime(G, varargin)

 [x,y] = GetClusterXYData(DATA, p)

 x = Rprime(r)

 Cut = PlotHistogram(DATA, E, varargin)

 ApplyLayout(DATA,varargin)

 [F, isnew] = SetFigure(lb, varargin)

 XYplot(a,b,tag,fcn)

 [x, xid] = GetValues(DATA, name, dim)

 PlotOneXY(DATA,names, varargin) 

 AddParameterMenu(F, callback, tag)

 CompareMean(a,b, p)

 DATA = LoadCellFile(DATA)

 [true, cellid] = isacell(DATA, ei, p)

 AddCellMenu(DATA)

 exitallv(src, evnt)

 vpts = SetVsamples(DATA, probe, np, nv)

 SetADC(a,b,fcn)

 ShowADCPos(src, data, type)

 CalcXcorr(a,b,fcn)

 FullVKeyPressed(src, ks)

 PCKeyPressed(src, ks)

 SelectTrial(src, b)

 c = TrialMarkChar(T)

 PlotOneTrial(DATA,id)

 idlist = SetTrialList(DATA)

 HistKeyPressed(src, ks)

 pos = RotateLine(pos, da)

 KeyPressed(src, ks)

 PlotGridSpikes(DATA, nspk, varargin)

 PlotQuickSpikes(DATA, nspk, varargin)

 [nt, spklst] = PlotTrialSpikes(DATA, nt, varargin)

 PlotCluster(a,b, mode)

 HistMenu(a,b, mode)

 ExptFigMenu(a,b, mode)

 PlotXcorr(a,b, pa, pb)

 Q = Cell2Cluster(cell,Clusters)

 SyncSpikes(DATA, cells)

 SpikeDraw(a,b, mode)

 xc = ShapeCorr(P,Q)

 cells = PlotAllXCorr(DATA, DataClusters, cells, varargin)

 ReplotXcorrs(a,b, type)

 PlotClusterXY(DATA, C)

 h = AddCellLabel(DATA,e,p)

 PlotAllProbes(DATA,type, varargin)

 c = MarkAxes(ax, mark)

 SummaryHit(a,b, p)

 PlotAllMeans(DATA)

 PlotMeanSpikes(C, p, cluster, varargin)

 bad = BadCluster(C)

 DATA = UseAllEvents(DATA)

 HitXYPlot(a,b, p)

 SpikeButtonPressed(a,b)

 HitImage(a,b,p)

 SpoolAllSpikes(DATA, varargin)

 PlotProbeSpikes(DATA, Vall, p, spklist,npts,offset)

 SpoolSpikes(DATA, varargin)

 csd = GetCSD(DATA, ndiff)

 ScrollV(src, evnt)

 ScrollSpikes(src, evnt)

 PlotFeatures(DATA, a, b, type, id, colors, C, varargin)        

 PlotPCs(pcs, a,b, type, id, colors, C,varargin)

 r = CalcRadius(E,xy)

 PlotVals(DATA, a,b, type, id, colors, varargin)

 DATA = ClassifyAll(DATA, force,varargin)

 need = NeedClusterData(Cluster, ci)

 D = CondenseCluster(C)

 DATA = IterateFit(DATA, niter)

 DATA = SetTemplateData(DATA,  cl, varargin)

 DATA = ClassifyAndFit(DATA)

 cluster = ClassifyFit(DATA, E, cnum)

 r = CalcClusterDistance(cluster, xy)

 X = GetDataStruct(DATA, f)

 [cl, cluster, xy] = ClassifySpikes(DATA, E, varargin)

 QuickSpks(DATA,nspk)

 cluster = PlotTriggerHist(DATA, cluster,varargin)

 chspk = UseProbeList(DATA, nprobes)

 MeanSpike = PlotMeanSpike(DATA, varargin)

 [Scores, T, details] = CalcScores(DATA, MeanSpike)

 Labels = PCLabels(DATA, usestd)

 Labels = TemplateLabels(DATA, usestd)

 [out, TemplateUsed, DprimeUsed] = TemplatePlot(DATA, varargin)

 [bs, as] = CalculateTemplateDips(DATA)

 DATA = ReplotPCs(DATA,E, varargin)

 clusterplot = GetClusterPlots(DATA,E, plots,pt)

 h = SetClusterIcon(DATA)

 clusterplot = GetClusterPlot(DATA,E,plots, pt)

 AddMarkToPCs(pos, space, plots, varargin)

 PlotVarE(DATA)

 DATA = RestrictTimeRange(DATA, t)

 DATA = ExcludeTrials(DATA, trials, add)

 SetCellCompare(a,b, cellid)

 AddCellMean(DATA, cellid)

 FinishSpikePlot(DATA)

 myAllV = GetAllV(DATA)

 ph = PlotSpikes(DATA,spkid, varargin)

 ShowFullV(src,b, fcn)

 [Vall, id] = PlotFullV(DATA, t, varargin)

 OldSetMenuCheck(F, tag, value)

 SetGUI(DATA)

 PlotTemplateScores(DATA, TemplateScores, probes)

 [dp, res] = MaxDprime(x, varargin)

  sgn = CheckSign(C, x, energy)

 [dip, details] = oldFindDip(values, energy, varargin)

 dp = CalcDprime(x, y)

 cname = ClusterFile(name,varargin)

 HistButtonPressed(src, data)

 HistButtonReleased(src, data)

 C = ClusterInfo(DATA)

 HistButtonDragged(src, data)

 h= oldDrawEllipse(E,varargin)

 h= oldDrawLine(E,varargin)

 h= DrawEllipse(E,varargin)

 h= DrawLine(E,varargin)

 [distance, cluster] = FindNearestCluster(DATA, pos)

 ButtonPressed(src, data)

 len = LineLength(l)

 distance = DistanceToEllipse(E, pos);

 in = InGraph(pt, ax)

 ButtonReleased(src, data)

 ScrollWheel(src, evnt)

 DATA = RotateCluster(DATA, angle)

 PCButtonDragged(src, data)

 pos = xyr2pos(xyr)

 C = GetClusterDef(cluster, n)

 DATA = SetEllipseDrawing(DATA, shape,varargin)

 MiscMenu(a, b, type)

 SetOption(a, b, type)

 PreferencesPopup(a,b, fcn)

 value = Text2Val(F, tag)

 PlotResult(a, b, type)

 HitTrial(data,b, cell)

 [Expt, plotres] = PlotExptCounts(DATA)

 PlotMenu(a, b, type, varargin)

 sdx = PlotXYSequence(DATA, probe, varargin)

 PlotISI(a, b, type)

 HitISI(a,b, t)

 [isis, trials, spkids] = CalcISI(Trials, varargin)

 [C, fits] = OptimizeEllipse(DATA)

 [SSD,dipr, details ] = EllipseDip(params, DATA, state)

 guess = MinimiseEllipseDip(C, DATA)

 [SSD, details ] = MinimiseEllipse(params, DATA, state)

 [SSD, details ] = MinimiseEllipseb(params, DATA, state)

 [C, fits] = OptimizeLine(DATA)

 Plot2GaussFit(params, DATA, state)

 [SSD, details ] = MinimiseLine(params, DATA, state)

 PlotGauss2Fit(fits);

 [dp, fits, details] = Fit2Gauss(C, r, DATA, varargin)

 DATA = LoadComments(DATA)

 AddComment(a,b,str)

 C= CondenseClusters(C, go, varargin)

 p = GetProbeFromName(name)

 ShowTaggedProbes(DATA)

 DATA = QuickAutoCut(a,b)

    end
end
