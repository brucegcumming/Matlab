% DEV
%
% Files
%   AddEllipse            - adds an ellipse to a figure interactively
%   AddLFPFT              - Adds Foutier Tranform to Trials in LFP
%   AddTimesToFullV       - Adds timestamp values for each Voltage sample to FullV
%   AddToList             - appends string 'name'  to text file 'list'
%   AlignMatrices         - [shift, xcs] = AlignMatrices(A,B, dim, varargin)
%   All2Expt              - Expt = All2Expt(AllExpt, cell, varargin)
%   AllExpt2mat           - [M, details] = AllExpt2Mat(E, ptype,varargin)
%   allSNRS               - find sessions with a spikes directory
%   AllVPcs               - AllVPcs(V, ...)  Extract and Classify Spikes from continuous voltage record
%   APlaySpkFile          - Lists Expts and Trials in .mat files from spike 2
%   AppDataPopup          - Builds a popup to edit elements of a structres
%   ApplyLayout           - DATA = ApplyLayout(DATA, varargin)
%   ArrayMovie            - ArrayMovie(Expt, ids, varargin)
%   BackupFile            - BackupFile(name)
%   BestAngles            - [theta, c, details] = BestAngles(x,y, varargin)
%   BlackRockPath         - Add Blackrock scripts to path
%   bondyadevpath         - 
%   BuildAllFullV         - Builds FullV Matlab files from component files exported by spike 2. 
%   BuildBWPsychRC        - [kernel, details] = BuildBWPsychRC(Expt, varargin)
%   BuildBWPsychRCa       - [kernel, details] = BuildBWPsychRC(Expt, varargin)
%   BuildFullVt           - builds vector of sample times from blkstart/blklen in a FullV struct
%   BuildGridFullV        - Called by ProcessGridFullV. Actuallly Builds the FullV files from .ns5
%   BuildGridIndex        - idx = BuildGridIndex(name, Expts, ...)
%   BuildGridIndexa       - Old way of building Grid index
%   BuildGridMeanSpike    - given a source probe, build spike triggered mean waveform on all probes
%   BuildMeanV            - Build MeanV file from individual Grid probe files.
%   BuildRFData           - Reads .ufl/rf.mat  files and builds table for PlotMatp
%   BuildRFFits           - find OP/PP expts in a dir and fit x,Y, positions
%   BuildSpkRecord        - Plots Block ON/Off and Trials for Expts
%   BWPsychPlot           - Plots a history of psychophysical performance given a range of dates,
%   CalcCellSize          - [result, details] = CalcCellSize(C, varargin) 
%   CalcCSD               - csd = CalcCSD(R, ...) calculates second derivative along dimension 1
%   CalcCSDPeaks          - 
%   CalcDprime            - 
%   CalcPostSaccGain      - 
%   CalcSubPlot           - x = function CalcSubPlot(nr,nc, probe, Array)
%   CalcTetrodeCorrs      - 
%   CalcVoffset           - V = CalcVoffset(AllV, chspk, overlap)
%   cellmatch             - id = cellmatch(s, C) returns and index of elements in C that match s
%   CellToMat             - CellToMat M = CellToMat(C, f,...
%   CellToStruct          - S = CellToStruct(C, varargin)
%   CheckClusterDir       - DATA = CheckClusterDir(name, varargin)  checks for any ClusterTimes that
%   CheckClusterLog       - 
%   CheckExceptions       - CheckExceptions(res, label)  runs through a result file
%   CheckExptRates        - 
%   CheckEyeData          - Expt = CheckEyeData(Expt) makes sure lengths of EyeData
%   CheckFullV            - CheckFullV(FullV, varargin)
%   CheckGridIndex        - 
%   CheckIdx              - 
%   CheckLamDatDir        - 
%   CheckLFP              - 
%   CheckLFPChoiceSign    - 
%   CheckList             - 
%   CheckPhysicalMemory   - find physical memory in Mb, using 'feature memstats'
%   CheckSpkMat           - CheckSpkMat(varargin)  checks contents of matlab variables from Spike3
%   ClearOSpike           - 
%   CloseAllHidden        - 
%   ClusterListIn230      - ClusterList('F:/Spike2/data/lem/M230','bycell',4,'landmarks',[5 18 46]);
%   ClusterLog            - ClusterLog(name,...) Show contents of ClusterLog
%   combine               - combine(file)
%   combine1              - combine(file)
%   combineA              - combine(file)
%   CombineCellsInList    - CombineCellsInList takes a CellList data file for combine and
%   CombineInit           - 
%   combineold            - combine(file)
%   CombineRFData         - rflist = CombineRFData(rflist, fits)
%   CombinerInit          - 
%   ComnineRFData         - rflist = CombineRFData(rflist, fits)
%   CompareAllClusters    - 
%   CompareAllExpt        - 
%   CompareClusterPair    - 
%   CompareClusters       - prints list of differnces between two Cluster structures
%   CompareExpts          - 
%   CompareFields         - [uniquea, uniqueb] CompareFields(a,b, .)
%   CompareGridIdx        - 
%   CompareIdx            - 
%   CompareNevNsx         - 
%   CompareSpikes         - If combine and AllVpcs are both up, compares spikes in the two.
%   CondenseClusters      - [CC, GM] = CondenseClusters(C, varargin)
%   ConvertClusterDetails - 
%   ConvertSpikeDir       - ConvertSpikeDir(path,type,varargin)
%   Copy of combine       - combine(file)
%   CorrMat               - CorrMat MATLAB code for CorrMat.fig
%   CorrMat1              - MATLAB code for CorrMat1.fig
%   cpfiles               - cpfiles (srcdir, tgtdir,...)
%   CreationDate          - Find the creation date of a file (E.g. smr file)
%   DEFINITIONS           - 
%   dir2name              - name = dir2name(path, type)
%   DrawClusterEllipse    - 
%   Expt2Name             - take an expt and return a name identifying the
%   ExptName              - 
%   ExptPlotXcorr         - 
%   ExptPsych             - function [pp, details] = ExptPsych(Expt, varargin)
%   ExptTrigLFP           - avg = ExptTrigLFP(Expt,LFP, varargin)
%   ExtractOnlineData     - 
%   FiguresToFront        - for all the tags conainted in the structure
%   FillLFPFields         - 
%   FilterFullV           - 
%   FindArtefacts         - 
%   FindDip               - Find a dip in a histogram of x;
%   FindEnclosingEllipse  - 
%   FindFrameTimes        - 
%   FitDriftMatrix        - 
%   FitSpike              - fit a simple descriptive spike shape such than when convolved with the 
%   FitSpikeShp           - fit a simple descriptive spike shape such than when convolved with the 
%   FixCluster            - FixCluster(C, varargin) modifids old clustr structs to make sure
%   FixClusterFile        - PlotClusters(dir, ......)
%   FixEd                 - 
%   FixEdepth             - 
%   FixExpt               - 
%   FixFullV              - 
%   FixLFP                - FixLPF removes the effect os a pike on the LFP
%   FixLFPFile            - LFP = FixLFPFile(name, varargin)
%   FixLFPMains           - LFP = FixLFPMains(LFP, tics)
%   FixLFPSpike           - Expt = FixLFPSpike(Expt)
%   FixLFPTrials          - Expt = FixLFPTrials(Expt)
%   FixMat                - 
%   FixSMR                - fixes errors in specifif matlab data files from Spike2
%   FixSpkDir             - 
%   FullV2LFP             - D = FullV2LFP(name,....)
%   FullVTimes2id         - id = VallTimes2id(Vall, t) 
%   GetArrayConfig        - ArrayConfig = GetArrayConfig(name,varargin)
%   GetExptNumber         - n = GetExptNumber extracts the number of an experiment from a filename
%   GetFilePath           - path GetFilePath(type,...)  returns path for datafiles, preferences, etc
%   gethostname           - 
%   GetMonkeyName         - [monk, monkyename dirsuffix] = GetMonkeyName(name)
%   GetPenInfo            - 
%   GetPenInfoStr         - txt = GetPenInfoStr(T, varargin)
%   GetPopStr             - getpopstr(x)  returns the curretly selected string value from a pop-menu
%   GetProbeFromName      - 
%   getSNR                - calculates the mean, std, and SNR for a 
%   GetUserName           - 
%   GMBoundary            - 
%   GMDip                 - [dip, details] = GMDip(xy, energy, varargin)
%   GMdprime              - calcualte drpime between two Gaussians in gmdistribution fit        
%   GMfit                 - [G, D, all] = GMfit(X, nd, nr, ..)
%   GoodCluster           - 
%   hostid                - 
%   HW5                   - 
%   InterpGauss1          - 
%   IsAllExpt             - 
%   KeepFigure            - Changes a figure Tag so that subsequent calls to GetFigure do not
%   LFPLatency            - [latency, details] = LFPLatency(LFP,times, varargin)
%   ListClusterBackup     - ListClusterBackup(prefix)
%   ListExpts             - [names, details] = ListExpts(E, ...)   lists names of Expts in  E
%   ListTree              - ListTree(Tree, pattern) Go through AGB disk Listings
%   LoadAllFullV          - 
%   LoadCellList          - 
%   LoadCluster           - LoadCluster(dirname, expts, ...)
%   LoadClusterFiles      - LoadClusterFiles(dir, ......)
%   LoadClusterInfo       - Expt = LoadClusterInfo(Expt, varargin)  Read additional
%   LoadClusters          - [AllClusters, AllFullVData] = LoadClusters(dirname, eid, varargin)
%   LoadFullV             - FullV = LoadFullV(name)
%   LoadGridLFP           - 
%   LoadLFP               - Expt = LoadLFP(Expt, ...)
%   loadmap               - map = loadmap(file)
%   MakeFullVInfo         - FullVData = MakeFullVInfo(Vall)
%   MakeProbeIndex        - probes = MakeProbeIndex(path, ...
%   MakeVisible           - 
%   ManData               - 
%   manSNRs               - find sessions with a spikes directory
%   MatchTimes            - 
%   MatrixPermute         - for square matrix M, re-oder rows and columnts according to order
%   MaxCohere             - 
%   MergeExcludedTrials   - 
%   MergeLFPExpts         - MergeLFPExpts(LFP)
%   mydialog              - %% Nargin Check %%%
%   MyDip                 - X = MyDip(X, varargin) Heuristic for finding best dip in a distribution
%   myscatter             - myscatter(x,y,symbol, ... makes a scatterplot where each datapoint has a callback
%   Name2Probe            - 
%   NewAllVPcs            - 
%   newPlotAllProbes      - PlotAllProbes(name, varargin) reads in Expt files for each probe, and
%   ns2nev                - 
%   OldBuildGridIndex     - idx = BuildGridIndex(name, Expts, ...)
%   OptimizeCluster       - 
%   Plot3Matrix           - 
%   PlotACResult          - PlotACResult(result) replots a result returned PlotExpt for AC expts, and
%   PlotAllCellFiles      - plot expts from an allexpt Structure
%   PlotAllExpts          - PlotAllExps plots results that combine multiple expts. 
%   PlotAllProbes         - PlotAllProbes(name, varargin) reads in Expt files for each probe, and
%   PlotAllXcorr          - 
%   PlotClusters          - PlotClusters(dir, ......)
%   PlotClusters1         - PlotClusters(dir, ......)
%   PlotComments          - Plotcomments(name...   Shows comments for an experimen/fiel/directory, and
%   PlotDataLine          - 
%   PlotExptFit           - Plots result from FitExpt
%   PlotFullV             - 
%   PlotGridFits          - 
%   PlotGui               - PlotGui(X,Y, Z, varargin)
%   PlotISI               - PlotISI(X, varargin)
%   PlotLFP               - 
%   PlotLFPprobes         - plot an LPF result stucture from PlotREvCOrAny 
%   PlotLFPpwr            - Plot responses generated by PlotMLFP;
%   PlotMap               - PlotMap('monkeyname')
%   PlotMeanSpike         - 
%   PlotMLFP              - [fpm, details] = PlotMLFP(Expt, ..)
%   PlotND                - PlotND(X,plots, varargin)
%   PlotOnline            - combine(file)
%   plotpen               - function plotpen(map, x, y)
%   PlotPlaidResult       - PlotPlaidResult(result) replots a result returned PlotExpt for Plaid expts, and
%   PlotProbes            - 
%   PlotQuickEMTrial      - 
%   PlotRates             - result = PlotRates(Expt,type, ...)
%   PlotRateSequence      - 
%   PlotRCLFPAll          - 
%   PlotSpikeC            - Build Structure with one large matrix of continuous voltages, from
%   PlotSpikes            - 
%   PlotSpikeShapes       - PlotSpikeShapes(name)
%   PlotSpikeTimes        - Plot Spike Times from Cluster lists made by AllVPCs
%   PlotSyncSpikes        - 
%   ProcessGridFullV      - DATA = ProcessGridFullV(name, varargin)
%   ProfileSummary        - ProfileSummary(P, varargin)  summarizes results from P = profile('info')
%   PulseTrigLFP          - 
%   RadialSum             - sor = RadialSum(kernel, fixr)
%   RC2mat                - x = RC2mat(rc, varargin)
%   ReadAutoLog           - 
%   ReadOnlineTxt         - read online data with just spike times in text format
%   ReadSerialFile        - [result, Expt, txt] = ReadSerialFile(name, 'readexpt')
%   ReadSerialLine        - 
%   ReadSpkFile           - 
%   ReadTemplateClusters  - ReadTemplateClustes(a)
%   ReadUfl               - res =ReadUfl(name)
%   ReCountExptSpikes     - 
%   RectAdd               - M = RectAdd(a,b,...) creates a matrix M from adding
%   RemoveFullVMains      - 
%   ReWriteClusters       - 
%   RewriteProbeSpikes    - 
%   rmfields              - rmfields(S, a,b,c....)    removes a list of fields from S
%   RunAllGridFiles       - RunAllGridFiles(file, .....
%   RunAllVPcs            - res = RunAllVPcs(dir, .....
%   RunPlotClusters       - 
%   SaveArrayConfig       - SaveArrayConfig(DATA) save array configuration
%   SaveLayout            - DATA = ApplyLayout(DATA, varargin)
%   ScaleLFP              - 
%   ScaleLFPRC            - 
%   SelectFits            - gid = SelectFits(fits,varargin)
%   SelectProbe           - selected = SelectProbe(probelist, selected)
%   SetClusters           - modifies ClusterTimesFiles on disk
%   SetFigPos             - SetFigPos(X, tag)
%   SetFigureName         - 
%   SetMenu               - 
%   SetPdir               - 
%   SetProbesToUse        - 
%   SetUIRecursive        - 
%   ShuffleTimes          - st = ShuffleTimes(t, Expt,varargin)
%   smooth2               - 
%   SpoolFullV            - AllVPcs(V, ...)  takes an MxN matrix (electrode by voltage) of continuous
%   SpoolSpikes           - GUI for playing back spikes from existing expts.
%   SummarizeCells        - 
%   TestEllipse           - 
%   TestEllipseDraw       - 
%   TestEvents            - 
%   TestFullV             - 
%   TestSpkFile           - [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%   TestSpkFileA          - [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%   TestSpkFileB          - [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%   TestStimOn            - TestStimOn
%   TestText              - 
%   TrackBlanks           - 
%   TrackLFPMains         - TrackLFPMains  adjusts an  LFP file for flucutating mains noise.
%   TwoProbes             - [aid, bid] = TwoProbes(sa,sb, varargin)
%   untitled2             - 
%   UpdateBabProbeList    - FInds probes markded as bad in Clusters, and Adds these to the ArrayConfig
%   VallTimes2id          - id = VallTimes2id(Vall, t) 
%   xcorrtimes            - [xc, details] = xcorrtimes(a, b, varargin)
%   XYClassify            - 
%   xyrotate              - XY = xyrotate(x,y,angle)
%   GetStructVal          - Get a value from a structure, checking fields exist.
