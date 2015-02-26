% MATLAB
%
% Files
%   ._iscluster             - 
%   acknowledge             - acknowledge(msg, windowlabel)  is just a wrapper for msgbox
%   AddEllipse              - adds an ellipse to a figure interactively
%   AddError                - DATA = AddError(DATA, varargin) prints an error and adds it to DATA struct
%   AddExpts                - Combines a list of Expts into one,
%   AddLFPFT                - Adds Foutier Tranform to Trials in LFP
%   AddRCFrameTimes         - Adds Start/End Times to an RC Expt
%   AddRPlot                - [rax,h] = AddRPlot(lax, varargin)
%   AddTimesToFullV         - Adds timestamp values for each Voltage sample to FullV
%   AddToList               - appends string 'name'  to text file 'list'
%   aliclicked              - buttonpress for alis figures
%   AlignMatrices           - [shift, xcs] = AlignMatrices(A,B, dim, varargin)
%   All2Expt                - Expt = All2Expt(AllExpt, cell, varargin)
%   AllExpt2mat             - [M, details] = AllExpt2Mat(E, ptype,varargin)
%   AllProbes               - Make sure lfpblank resp is included when doing all files.
%   allSNRS                 - find sessions with a spikes directory
%   AllVPcs                 - AllVPcs(V, ...)  Extract and Classify Spikes from continuous voltage record
%   angle3d                 - x = angle3d(a,b)
%   anova1u                 - ANOVA test whether it is tuned to disparity
%   APlaySpkFile            - Lists Expts and Trials in .mat files from spike 2
%   AppDataPopup            - Builds a popup to edit elements of a structres
%   ApplyLayout             - DATA = ApplyLayout(DATA, varargin)
%   ArrayDistance           - 
%   ArrayMovie              - ArrayMovie(Expt, ids, varargin)
%   autocorrelate           - [ACF, ACFeach, nspikes] = autocorrelate(Trials, latency, duration)
%   autocorrelate_bins      - [ACF ACFeach] = autocorrelate(Trials, latency, binw)
%   BackupFile              - BackupFile(name)
%   BenchMat2               - matlab benchmark meant to resemble LSR/Vision work
%   BestAngles              - [theta, c, details] = BestAngles(x,y, varargin)
%   BestF1F0                - 
%   bgcload                 - load matlab data, and include pathname
%   BinocCommand            - BinocCommand(str) Sends a single commnad to binoclean
%   binocmovie              - make avi file out of image frames from Binoc
%   BlackRockPath           - Add Blackrock scripts to path
%   Bresample               - Resample(values,nrep,nout)
%   Bruce                   - Functions written by Bruce for Expt .mat files:
%   BuildAllFullV           - Builds FullV Matlab files from component files exported by spike 2. 
%   BuildBWPsychRC          - [kernel, details] = BuildBWPsychRC(Expt, varargin)
%   BuildBWPsychRCa         - [kernel, details] = BuildBWPsychRC(Expt, varargin)
%   BuildExptComments       - BuildExptComments(path, varargin) Extracts comment lines from
%   BuildFileName           - str = BuildFileName(name, type, ...) returns string naming file of type
%   BuildFullVt             - builds vector of sample times from blkstart/blklen in a FullV struct
%   BuildGridFullV          - Called by ProcessGridFullV. Actuallly Builds the FullV files from .ns5
%   BuildGridIndex          - idx = BuildGridIndex(name, Expts, ...)
%   BuildGridIndexa         - Old way of building Grid index
%   BuildGridMeanSpike      - given a source probe, build spike triggered mean waveform on all probes
%   BuildMeanV              - Build MeanV file from individual Grid probe files.
%   BuildPath               - path = BuildPath(name, varargin)
%   BuildRFData             - Reads .ufl/rf.mat  files and builds table for PlotMatp
%   BuildRFFits             - find OP/PP expts in a dir and fit X,Y, positions
%   BuildSpkRecord          - Plots Block ON/Off and Trials for Expts
%   BWPsychPlot             - Plots a history of psychophysical performance given a range of dates,
%   CalcCellSize            - Reads Clusters (AllVPcs) and calculates size, and vertical spread
%   CalcClust               - 
%   CalcClusterVars         - DATA = CalcClusterVars(DATA, ispk, varargin)
%   CalcCP                  - function [cp, details] = CalcCP(acounts, bcounts, varargin)
%   CalcCSD                 - csd = CalcCSD(R, ...) calculates second derivative along dimension 1
%   CalcCSDPeaks            - res = CalcCSDPeaks(csd, varargin) Find Peak time in csd response
%   CalcDDI                 - ddi = CalcDDI(Expt, duration)
%   CalcDprime              - dp = CalcDprime(x, y) Calc dprime for two distributions x and y
%   CalcEfficacies          - 
%   CalcEfficacy            - [e, nspk] = CalcEfficacy(x,y,varargin)
%   CalcFAmiss              - 
%   CalcFixDisp             - 
%   CalcGrandCP             - 
%   CalcISI                 - [isis, trials, spkids] = CalcISI(Trials, varargin)
%   CalcIsolation           - [x, details] = CalcIsolation(pts, idlist, clnum...   calculate islolation metrics from data points
%   CalcLFPPulse            - [lfp,n] = CalcLFPPulse(Expt, All, ...)
%   CalcPk                  - Calculate a psychophysical kernel from revcor Expt data
%   CalcPostSaccGain        - Analyses OPRC data perisaccadically
%   CalcRebound             - 
%   calcsacsdf              - calcsdf takes a vector list of spike times, spikes
%   calcsdf                 - calcsdf takes a vector list of spike times, spikes
%   CalcSdfs                - 
%   CalcSubPlot             - x = function CalcSubPlot(nr,nc, probe, Array)
%   CalcTetrodeCorrs        - 
%   CalcVoffset             - V = CalcVoffset(AllV, chspk, overlap)
%   calibrateICCStereo      - oldICCFilePath the path to the icc profile used to do the measurements
%   CatStruct               - Both = CatStruct(To, From) Add struct elements even if fields don't match
%   CellFields              - CellFields(C, field, varargin) finds elements of C contianing field f
%   CellIsEmpty             - true CellIsEmpty(S, field, n) Checks struct S to see if cell array
%   cellmatch               - id = cellmatch(s, C) returns and index of elements in C that match s
%   cellmember              - cellmember(a,b) ismember for cell arrays of doubles.
%   cellregexp              - cellstrfind(str, pat) calls strfind on each cell in a cell array str, to
%   cellstrfind             - cellstrfind(str, pat) calls strfind on each cell in a cell array str, to
%   CellToMat               - CellToMat M = CellToMat(C, f,...
%   CellToStruct            - S = CellToStruct(C, varargin)
%   cf                      - 
%   CheckAllExpts           - 
%   CheckAutoCuts           - 
%   CheckClusterDir         - DATA = CheckClusterDir(name, varargin)  checks for any ClusterTimes that
%   CheckClusterLog         - 
%   checkclustertimes       - 
%   CheckExceptions         - CheckExceptions(res, label)  runs through a result file
%   CheckExpt               - isi = CheckExpt(Expt)
%   CheckExptClusters       - CheckExptClusters(D, varargin) checks for misisng trials in Clusters
%   CheckExptRates          - [err, counts] = CheckExptRates(Expt, varargin)
%   CheckExpts              - [good, res] CheckExpts(Expts, varargin)
%   CheckEyeData            - Expt = CheckEyeData(Expt) makes sure lengths of EyeData
%   CheckFileUpdate         - CheckFileUpdate(src, tgt, chkmode) Copy new files
%   CheckForPCA             - 
%   CheckFrameDiffs         - Gnve a set of frametime differences, find skips
%   CheckFrameTimes         - Go through Trials from bnc file checking for dropped frames
%   CheckFullV              - CheckFullV(FullV, varargin)
%   CheckGridIndex          - 
%   CheckIdx                - 
%   CheckLamDatDir          - 
%   CheckLFP                - 
%   CheckLFPChoiceSign      - 
%   CheckList               - 
%   CheckMatlab             - CheckMatlab(name) Count number of lines in functions for a matlab source file 
%   CheckNameBug            - name = CheckNameBug(name)
%   CheckPath               - out = CheckPath(path) tries to ensure that paths for data files are corret
%   CheckPhysicalMemory     - find physical memory in Mb, using 'feature memstats'
%   CheckPtSize             - 
%   CheckRespDir            - returns the sign convention of respdir, based on the 
%   CheckSaccades           - 
%   CheckSerialFile         - Check for errors in serial file, like failures to set dx in manual expts
%   CheckSpike              - 
%   CheckSpikeCh            - 
%   CheckSpikeChForValues   - CheckSpikeChForValues
%   CheckSpikePair          - 
%   CheckSpk2               - 
%   CheckSpkFile            - Check spk files are up to date for a cluster
%   CheckSpkMat             - CheckSpkMat(varargin)  checks contents of matlab variables from Spike3
%   CheckStimDur            - [durs, details] = CheckStimDur(name,varargin) Read online serial file and checks stimulus durations 
%   checktimeranges         - given a directory, such as 'F:/Spike2/data/lem/M195'
%   CheckTrees              - [missing, details] = CheckTrees(src,tgt) check that files on src exist somewhere on tgt
%   circgauss               - circgauss(params,x) evalutates a circular Gaussian function of an input x, in radians.
%   circle                  - h = circle(r, c, args) plots a circle, with args{:} passed on to plot()
%   CleanSpikes             - CleanSpikes(ch, ...
%   ClearOSpike             - 
%   ClearPlot               - ClearPlot() Clears a plot - deleting all axes 
%   CloseAllHidden          - 
%   CloseChildren           - CloseChildren(parent, closes figures associated with parent
%   CloseTag                - 
%   ClusterDef              - returns list of defined clusters in C
%   ClusterDistance         - convert Cluster.xy to radius/angle
%   Clustering              - List of Matlab functions/structure for Clustering
%   ClusterIsSet            - --- help for AllV/ClusterIsSet ---
%   ClusterListIn230        - ClusterList('F:/Spike2/data/lem/M230','bycell',4,'landmarks',[5 18 46]);
%   ClusterLog              - ClusterLog(name,...) Show contents of ClusterLog
%   ClusterParam            - Cluster Structures created from FullV Files
%   ClusterStrings          - 
%   CmpSMRs                 - CmpSMRs(monkey) finds files for monkey on C:\Spike6\data and on Z:\smr
%   Code2String             - Code2String(code, type) return error string from codes
%   ColorZoom               - ColorZoom(zoom)  Zoom in/out a color axis. >1 = produce larger color range
%   combine                 - cmb.combine(file)
%   combine1                - combine(file)
%   combineA                - combine(file)
%   CombineAllExpts         - 
%   CombineAutoCut          - 
%   CombineCalcDensity      - [x,y,z] = CalcDensity(DATA, expspks, mode)
%   CombineCellsInList      - CombineCellsInList takes a CellList data file for combine and
%   CombineExpts            - 
%   CombineInit             - 
%   combineold              - combine(file)
%   CombineRFData           - rflist = CombineRFData(rflist, fits)
%   CombinerInit            - 
%   ComnineRFData           - rflist = CombineRFData(rflist, fits)
%   CompareAllClusters      - 
%   CompareAllExpt          - 
%   CompareClusterPair      - res = CompareClusterPair(A,B, varargin) Compares mean spike shape of two Clusters
%   CompareClusters         - prints list of differnces between two Cluster structures
%   comparedates            - compares files in two directories.  For files with the same name, it
%   CompareDetails          - CompareDetails(dpath, varargin) Check folder dpath for erros in
%   CompareExpts            - CompareExpts(Ea,Eb, varargin) Compares Expt struct and reports differences
%   CompareFields           - [uniquea, uniqueb] CompareFields(a,b, .)
%   CompareGridIdx          - 
%   CompareIdx              - 
%   CompareLists            - Compare contents of two lists of files
%   CompareNevNsx           - 
%   CompareSpikes           - If combine and AllVpcs are both up, compares spikes in the two.
%   CompareSpikeWaves       - CompareSpikeWaves(P,Q,aid,bid) Compare a set of spikes 
%   compdirectory           - GUI tool for directory comparison.
%   CondenseClusters        - [CC, GM] = CondenseClusters(C, varargin)
%   CondenseDetails         - 
%   CondenseRC              - CondenseRC(Expt, code)
%   confirm                 - confirm(str) calles questdlg, returns 0 or 1
%   confrim                 - confirm(str) calles questdlg, returns 0 or 1
%   ConsecUpper             - returns the last occurence of 2 or more consecutive upper-case
%   ConvertClusterDetails   - 
%   ConvertSpikeDir         - ConvertSpikeDir(path,type,varargin)
%   ConvertTime             - truet = ConvertTime(Expt, t) Converts timestamps into matlab datenums
%   Copy of combine         - combine(file)
%   CopyFields              - to = CopyFields(to, from)
%   CopySFields             - to = CopySFields(to, element, from) fopy fileds to structure element
%   CopyUiProperties        - args = CopyUiProperties(it)
%   corr_counts             - [counts, rnd]=corr_counts(r,ntrials,poolsize,varargin)
%   CorrectAuto             - 
%   CorrMat                 - CorrMat MATLAB code for CorrMat.fig
%   CorrMat1                - MATLAB code for CorrMat1.fig
%   CorrMat8G               - neurons in each pool
%   CountClusterList        - 
%   CountClusters           - 
%   Counts                  - [counts, vals] = Counts(x) 
%   cpfiles                 - cpfiles (srcdir, tgtdir,...)
%   cprintf                 - displays styled formatted text in the Command Window
%   CreationDate            - Find the creation date of a file (E.g. smr file)
%   cv2sd                   - returns SD of a cicular Gaussian that produces a give Circular Variance
%   DataThief               - DataThief(im,....)
%   DEFINITIONS             - 
%   DensityPlot             - handles = DensityPlot(x,y,varargin)
%   DeTrendExpt             - 
%   devpath                 - 
%   dir2name                - name = dir2name(path, type) superceded by BuildPath
%   DirCompare              - 
%   DiskSummary             - result = DiskSummary(path)
%   DistanceMeasure         - 
%   DrawClusterEllipse      - 
%   DrawEllipse             - DrawEllipse(xyr, varargin) draws ellipse
%   EmAvg                   - 
%   expand                  - 
%   ExportPlots             - ExportPlots(a,b, caller)
%   Expt2Blocks             - Expts = Expt2Blocks(E, varargin) Convert a combined Expt into a cell array of individual Blocks
%   Expt2Images             - Rebuild RDS patterns used in an Experiment
%   Expt2Name               - take an expt and return a name identifying the
%   ExptCellQuality         - 
%   ExptCorr                - [exptxc, details] = ExptCorr(exa,exb,conditions,varargin)
%   ExptListCorrs           - ExptListCorrs(list....
%   ExptName                - 
%   ExptPlotXcorr           - 
%   ExptPsych               - function [pp, details] = ExptPsych(Expt, varargin)
%   ExptSpikeListAll        - 
%   ExptString              - str = ExptString(E, s) Return a string descibing property s of Expt
%   ExptStruct              - Expt(Expt) - help on the Expt Structure
%   ExptTrigLFP             - avg = ExptTrigLFP(Expt,LFP, varargin)
%   ExtractOnlineData       - 
%   Extrema2d               - [maxs, mins] = Extrema2d(Z, ...)
%   EyeCal                  - [gain, pos, stds] = EyeCal(Expt, varargin)
%   EyePosAcov              - [ACF, ACFeach, avg] = EyePosAcov(Trials, ...)
%   famp                    - [a, c] = famp(x,y,f)
%   fgf                     - bring current figure to the front
%   figify                  - figify(fighandle, axhandle)  makes a matlab figure better for importing to PowerPoint, 
%   FigureForBruce          - 
%   FigureProp              - return figure number indepdent of matlab version
%   FiguresToFront          - for all the tags conainted in the structure
%   FillLFPFields           - 
%   fillpmesh               - [X , Y, Z] = fillpmesh(x,y,z)
%   FillTrials              - Expt = FillTrials(Expt, code, ...)
%   FilterFullV             - 
%   filterim                - filterim( sf, or, rsd, varargin)
%   FilterLFP               - [lfp, ft] = FilterLFP(lfp, samplerate)
%   FilterLFPTrials         - 
%   FindArtefacts           - 
%   FindChild               - calls find obj for a window and any paired windows
%   FindDip                 - Find a dip in a histogram of x;
%   FindDTFiles             - 
%   FindEllipse             - [x, score] = FindEllipse(xy,idlist, varargin) find the ellipse that encloses points
%   FindEnclosingEllipse    - 
%   FindExpt                - FindExpt(idxfile,'name1',name2,....)
%   FindFig                 - 
%   FindFiles               - FindFiles(monk) finds all usefule data directorys for an animal
%   FindFrameTimes          - 
%   FindLayerIV             - 
%   findmid                 - find center of a kernel
%   FindMissingTrials       - [badtrial,badid] = FindMissingTrials(Expt, t, varargin)
%   FindPenLog              - name = FindPenLog(monkey, pen, ...
%   FindSpikes              - [ispk, spktimes, codes] = FindSpikes(DATA, times, probe, range)
%   FindTrial               - tid = FindTrial(Expt, t, varargin) find Trial that contains time t
%   FitACGabors             - function para = Fit_Gabor(x,y,init)
%   FitDriftMatrix          - 
%   FitExp                  - Fit = FitExp(x,y,....)  fits an exponential y = exp(-x/a)
%   FitExpt                 - FitExpt(Expt....)
%   FitGabor                - =================================================================================================
%   FitGauss                - [FitGauss,fval,exitflag] = FitGauss(x,y, varargin)
%   FitOriTuning            - [FitGauss,fval,exitflag] = FitOriTuning(x,y, varargin)
%   FitPhaseGauss           - FitPhaseGuass(x,y) fits the sum of a Gaussian and a Cumulative Gaussian
%   FitPhaseGaussAll        - FitPhaseGuass(x,y) fits the sum of a Gaussian and a Cumulative Gaussian
%   fitpsf                  - function results = fitpsf(probitin, varargin)
%   FitSine                 - Fit = FitSine(x, y, varargin)
%   FitSpike                - fit a simple descriptive spike shape such than when convolved with the 
%   FitSpikeShp             - fit a simple descriptive spike shape such than when convolved with the 
%   FitSqrtGabor            - =================================================================================================
%   FitTopography           - FitTopograpy(x,y,px,py) Fit RF topography for penetraions
%   fittype                 - Fittype for curve and surface fitting
%   FixAllExpt              - FixAllExpt(A) Make shure AllExpt Stuct is consistent
%   FixCluster              - FixCluster(C, varargin) modifids old clustr structs to make sure
%   FixClusterFile          - FixClusterFile(dir, ......) Read Clusters, check, modify if necessary
%   FixClusterTimes         - fix mismatch between Clusters{c}.times and ClusterDetails{c}
%   FixEd                   - 
%   FixEdepth               - 
%   FixExpt                 - Expt = FixExpt(Expt,type)  Fixes misleading fields in Expt Structs
%   FixFullV                - 
%   FixLFP                  - FixLPF removes the effect os a spike on the LFP
%   FixLFPFile              - LFP = FixLFPFile(name, varargin)
%   FixLFPMains             - LFP = FixLFPMains(LFP, tics)
%   FixLFPSpike             - Expt = FixLFPSpike(Expt)
%   FixLFPTrials            - Expt = FixLFPTrials(Expt)
%   FixMat                  - 
%   FixSerialExpt           - Expt = FixSerialExpt(Expt)
%   FixSerialFile           - Use a serial File to fix missing info for a file
%   FixSMR                  - fixes errors in specifif matlab data files from Spike2
%   FixSmrMat               - FixSmrMat(name, ..) Fix errors in a .mat file made by Spike2.
%   FixSpkBlk               - FixSpkBlk(namestr, fixtime)
%   FixSpkDir               - 
%   FixSpkMains             - pavg = FixSpkMains(spk, tics, ...)
%   FullV2LFP               - D = FullV2LFP(name,....)
%   FullVTimes2id           - id = VallTimes2id(Vall, t) 
%   Gabor                   - Gabor(params, ...) generates a Gabor function, in one or 2 dimensions
%   gabor2d                 - 
%   Gauss                   - gauss(params,x) evalutates a Gaussian function of an input x
%   gauss2d                 - [X,Y,Z] = gauss2d(params,x) evalutates a Gaussian function of an input x
%   GaussStr                - 
%   GetArrayConfig          - ArrayConfig = GetArrayConfig(name,varargin)
%   GetCellNumber           - Get a cell number from a string,Expt or Cluster
%   GetCellValue            - GetCellValue(C, n, field) returns C{n}.field IF
%   GetDataFromFig          - if a subrouting has been called from a uicontrol, this find
%   GetElectrodeName        - GetElectrodeName(Expt, varargin) find string describing electrode in Expt
%   GetEval                 - [value, list] =  GetEval(Expt,type,flag)
%   GetExptNumber           - n = GetExptNumber extracts the number of an experiment from a filename
%   GetExptSpikeTimes       - 
%   GetFigPos               - GetFigPos(parent, varargin) find current locations of all child figures
%   GetFigure               - fign =Getfigure(tag,...)   is a convenience routine that finds a figure
%   GetFileName             - name = GetFileName(fpath) Retuns name part of a string (after directory)
%   GetFilePath             - path GetFilePath(type,...)  returns path for datafiles, preferences, etc
%   GetFrameReversal        - GetFrameReversal(mlfp, lfptimes, probes)
%   GetFullVName            - 
%   gethostname             - 
%   GetLFPSamples           - LFP = GetLFPSamples(Expt, times) Extract a matrix of LFP values from
%   GetMonkeyName           - [monk, monkyename dirsuffix, xname] = GetMonkeyName(name) - extract monkey, folder from name
%   GetName                 - Getname(X) returns the cell name associated with string/structre X
%   GetNames                - GetNames(file, exptid, suffix, suffix,.....
%   GetOptionalArg          - [j, value] = GetOptionalArg(args,j, default) get number from next arg in list,
%   GetParentFigure         - 
%   GetPenData              - 
%   GetPenInfo              - 
%   GetPenInfoStr           - txt = GetPenInfoStr(T, varargin)
%   GetPopStr               - getpopstr(x)  returns the curretly selected string value from a pop-menu
%   GetProbe                - pid = GetProbe(DATA, eid, probe)
%   GetProbeFromName        - 
%   GetProbeSep             - 
%   getSNR                  - calculates the mean, std, and SNR for a 
%   GetSpikeVals            - [x,DATA] = GetSpikeVals(DATA, ispk, values, dVdt, type, recalc, pcs)
%   GetSpkTimes             - extract spike times from and expt struct.
%   getstr2num              - val = getstr2num(tag) 
%   GetString               - f = GetString(tag, parent, callback)
%   GetStructVal            - Get a value from a structure, checking fields exist.
%   GetTextVal              - GetTextVal(args) used findobj(args) to locat a text GUI field,
%   GetTF                   - 
%   GetTFNames              - 
%   GetTriggerSet           - List  which cluster belongs to which triggerset
%   GetUserName             - 
%   git                     - A thin MATLAB wrapper for Git.
%   GMBoundary              - 
%   GMDip                   - [dip, details] = GMDip(xy, energy, varargin)
%   GMdprime                - calcualte drpime between two Gaussians in gmdistribution fit        
%   gmean                   - 
%   GMfit                   - [G, D, all] = GMfit(X, nd, nr, ..)
%   GoodCluster             - 
%   GrandCP                 - 
%   GroupByDistance         - Given a matix of mutual overlaps, break into groups
%   GroupByOverlap          - Given a matix of mutual overlaps, break into groups
%   HartigansDipSigniftest  - function		[dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)
%   HartigansDipTest        - function	[dip,xl,xu, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf)
%   here                    - Useful functions in My Matlab path
%   hostid                  - 
%   HW5                     - 
%   HW5dim                  - Problem 3 plotting Ds rather than gammas and lambdas
%   IDString                - 
%   InteractiveEllipse      - Draws an ellipse interactively
%   Interpf                 - Z = interpf(x, y, z, X, Y, type, smooth)
%   InterpGauss1            - 
%   InterpretLine           - InterpretLine(s, ....)
%   IsAllExpt               - 
%   iscluster               - res = iscluster(C) determine if C is a cluster strcutre
%   IsDataLine              - 
%   isfigure                - function yesno = isfigure(h)
%   jcontrolDemo            - jcontrolDemo demo file to illustrate some uses of jcontrol objects
%   KeepFigure              - KeepFigure(figure,  ...) Changes a figure Tag so that subsequent calls to GetFigure do not
%   LastConsecUpper         - returns the last occurence of 2 or more consecutive upper-case
%   LFPCohere               - 
%   LFPGains                - 
%   LFPLatency              - [latency, details] = LFPLatency(LFP,times, varargin)
%   listadd                 - listadd(varagin)
%   ListClusterBackup       - ListClusterBackup(prefix)
%   ListExpts               - [names, details] = ListExpts(E, ...)   lists names of Expts in  E
%   listop                  - listop(list, 'copy', ....)
%   ListTree                - ListTree(Tree, pattern) Go through AGB disk Listings
%   LoadAllClusterFiles     - 
%   LoadAllFullV            - 
%   LoadBwLFP               - 
%   LoadCellList            - 
%   LoadCluster             - LoadCluster(dirname, expts, ...)
%   LoadClusterDetails      - [ClusterDetails, details] = LoadClusterDetails(name)
%   LoadClusterFiles        - LoadClusterFiles(dir, ......)
%   LoadClusterInfo         - Expt = LoadClusterInfo(Expt, varargin)  Read additional
%   LoadClusters            - [AllClusters, AllFullVData] = LoadClusters(dirname, eid, varargin)
%   LoadComments            - 
%   LoadEmData              - [Expt, details] = LoadEmData(Expt) finds the matching .em matlab file, read in
%   LoadErrors              - Load error files from disk
%   LoadExpt                - Expt = LoadExpt(name,...
%   LoadFullV               - FullV = LoadFullV(name)
%   LoadGridLFP             - Expt = LoadGridLFP(Expt, varargin) loads LFP data into Expt with Utah Data
%   LoadLFP                 - Expt = LoadLFP(Expt, ...)
%   loadmap                 - map = loadmap(file)
%   LoadRefClusters         - 
%   LoadSpike2LFP           - [Expt, details] = LoadSpike2LFP(Expt, varargin)
%   MakeAllExpts            - MakeAllExpts(name, varargin) reads in Expt files for each probe, and
%   makedefaultkernels      - 
%   MakeDummyTrials         - Trials = MakeDummyTrials(...)
%   MakeFullVInfo           - FullVData = MakeFullVInfo(Vall)
%   makenewORBWRCkernels    - Searches through psych folders @
%   makenewORBWRCkernels___ - Searches through psych folders @
%   MakeProbeIndex          - probes = MakeProbeIndex(path, ...
%   MakeSpike2Paths         - Given pathname for a new Spike2 Data file, make sure necessary folders exist
%   MakeVisible             - 
%   ManData                 - 
%   manSNRs                 - find sessions with a spikes directory
%   MapDefs                 - 
%   MarkSlope               - 
%   MatchInd                - ymatch = MatchInd(x,y)
%   MatchTimes              - ids = MatchTimes(ta, tb, tw) find times in tb that match ta within +- dw
%   MatrixCorrCounts        - [counts, rnd]=corr_counts(r,ntrials,poolsize,varargin)
%   MatrixPermute           - for square matrix M, re-oder rows and columnts according to order
%   MaxCohere               - 
%   MeanVector              - r = MeanVector(lens, angles, ...)
%   MergeCluster            - takes a matrix of exptno/probe row vectors and copies the relevant
%   MergeExcludedTrials     - 
%   MergeLFPExpts           - MergeLFPExpts(LFP)
%   min2sec                 - 
%   minmax                  - 
%   mkpath                  - mkpath(path)
%   mksacsdf                - [sdf, nsac, nspikes, saccades] = mksacsdf(Trials, smooth, times, type, varargin)
%   mode                    - deprecated. Use prctile(x,50);
%   mssg                    - 
%   mupdate                 - update network/local copy of file from local/network copy
%   mycolors                - returns 126 color triplets, designed to be useful for plotting. The
%   mycprintf               - wrapper for cprintf that keeps currentfigure
%   mydialog                - %% Nargin Check %%%
%   MyDip                   - X = MyDip(X, varargin) Heuristic for finding best dip in a distribution
%   mydir                   - d = mydir(path)
%   myfittype               - FITTYPE   Fittype for curve and surface fitting
%   mygetCurrentTask        - like current task, but returns valid structure when no pool open.
%   myhandle                - myhandle(a) Checks that a is a graphics handle AND is not zero
%   myhostid                - 
%   myjcontrol              - 
%   mylegend                - mylegend(h, labels, varargin)
%   myNormalize             - Y = myNormalize(X,...)
%   myPCA                   - Do quick PCA on matrix, plot clustersing;
%   myplotyy                - 
%   myround                 - allows you to specify the number of decimal places to round to
%   myscatter               - myscatter(x,y,symbol, ... makes a scatterplot where each datapoint has a callback
%   mysubplot               - mysubplot(nr,nc,n, varargin)
%   mytoc                   - took mytoc(start)
%   name2path               - [path, dirpath] = name2path(name,...) Generate full path from name
%   Name2Probe              - 
%   NetFilename             - 
%   NewAllVPcs              - 
%   NewClusterStruct        - 
%   NewFigure               - 
%   newPlotAllProbes        - PlotAllProbes(name, varargin) reads in Expt files for each probe, and
%   NNdistance              - estimate cluster separation based on distance to neighbors
%   ns2nev                  - 
%   NsineRC                 - NsineRC(Expt, All, delays, ids, ....)
%   Nsinesdf                - res = Nsinesdf(Expt, im, varargin)
%   Nsubplots               - [nrow, ncol] = Nsubplots(nplots, varargin)
%   NumberFromString        - find pattern in string, and get following number
%   NumberToExptName        - 
%   ObjectVersion           - 
%   OldBuildGridIndex       - idx = BuildGridIndex(name, Expts, ...)
%   oldcombine              - combine(file)
%   op2xy                   - op2xy(op, angle,..)
%   OptimizeCluster         - 
%   OriginalAllVPcs         - AllVPcs(V, ...)  takes an MxN matrix (electrode by voltage) of continuous
%   OriginalPlotClusters    - PlotClusters(dir, ......)
%   ParseExptComment        - Expt = ParseExptComment(Expt,s)
%   parzenSurf              - parzenSurface -------------------------------------
%   parzenSurf2d            - parzenSurface -------------------------------------
%   path2name               - path2name(str...) take full path and return just distinctive name path
%   PaychMon                - 
%   PlaySpikes              - 
%   PlaySpkFile             - [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%   Plot3Matrix             - 
%   PlotACResult            - PlotACResult(result) replots a result returned PlotExpt for AC expts, and
%   PlotAllCellFiles        - plot expts from an allexpt Structure
%   PlotAllExpt             - plots data from an AllExpt Struct
%   PlotAllExpts            - PlotAllExps plots results that combine multiple expts. 
%   PlotAllProbes           - PlotAllProbes(name, varargin) reads in Expt files for each probe, and
%   PlotAllXcorr            - 
%   PlotCellIm              - PlotCellIm(CellList, CellDetails, nclusters, varargin)
%   PlotClusters            - PlotClusters(dir, ......)
%   PlotClusters1           - PlotClusters(dir, ......)
%   PlotComments            - Plotcomments(name...   Shows comments for an experimen/fiel/directory, and
%   PlotDataLine            - 
%   PlotEm                  - PlotEm(Expt, trials, ......)
%   PlotErrors              - PlotErrors(name...   Show errors and warnings for and expt directory
%   PlotEvents              - 
%   PlotExpt                - result = PlotExpt(name,'varargin')
%   PlotExptEM              - PlotExptEM plots eye movement data (and loads the data into Expt if
%   PlotExptFit             - Plots result from FitExpt
%   PlotExpts               - PlotExpts(dir, ...)
%   PlotExptSaccades        - 
%   PlotExptSeq             - Plot sequnece of events in an Expt. Default is to plot Stimon/off times
%   PlotExptSpikes          - PlotExptSpikes(Expt, varargin)
%   PlotExptsSummary        - PlotExptsSummary(Expts) Summarizes expts in a cell array or PsychFile 
%   PlotFiles               - runonline(list)
%   PlotFullV               - 
%   PlotGridFits            - 
%   PlotGui                 - PlotGui(X,Y, Z, varargin)
%   PlotISI                 - PlotISI(X, varargin)
%   PlotLFP                 - 
%   PlotLFPprobes           - plot an LPF result stucture from PlotREvCOrAny 
%   PlotLFPpwr              - Plot responses generated by PlotMLFP;
%   PlotMap                 - PlotMap('monkeyname')
%   PlotMeanSpike           - 
%   PlotMLFP                - [fpm, details] = PlotMLFP(Expt, ..)
%   PlotND                  - PlotND(X,plots, varargin)
%   plotonegrid             - theplot controls plot type. 
%   PlotOnePen              - PlotOnePen(penlog,varargin) plots data in a penetration log named in
%   PlotOnline              - combine(file)
%   PlotOnlineExpt          - plots online data delivered from spike2
%   PlotOXM                 - 
%   plotpen                 - function plotpen(map, x, y)
%   PlotPGM                 - 
%   PlotPlaidResult         - PlotPlaidResult(result) replots a result returned PlotExpt for Plaid expts, and
%   PlotProbes              - 
%   plotpsych               - 
%   PlotQuickEMTrial        - 
%   PlotRates               - result = PlotRates(Expt,type, ...)
%   PlotRateSequence        - PlotRateSequence(Expts, ...) plots the rates for each trial in a set of expts
%   PlotRC                  - PlotRC takes a result file from PlotRevCorAny, and replots the data
%   PlotRCLFP               - Plot rc struct returned by PlotRevCorAny, for LFP
%   PlotRCLFPAll            - 
%   PlotResult              - PlotResult(result) replots a result returned by PlotExpt()
%   PlotRevCor              - Add uncorr to plot
%   PlotRevCorA             - result = PlotRevCorAny(Expt....)
%   PlotRevCorAny           - result = PlotRevCorAny(Expt....)
%   PlotRFFits              - plots RF fits made with BuildRFFits
%   plotsacs                - Add option to show stimulus box.circle
%   PlotSdfs                - 
%   plotsdpsych             - 
%   PlotSMREvents           - 
%   plotspike               - 
%   PlotSpikeC              - Build Structure with one large matrix of continuous voltages, from
%   PlotSpikeFile           - PlotSpikeFile(file, varargin)
%   PlotSpikes              - 
%   PlotSpikeShapes         - plots mean spike shape for each Expt/Probe in psuedocolor.
%   PlotSpikeTimes          - Plot Spike Times from Cluster lists made by AllVPCs
%   PlotSpin                - 
%   PlotSyncSpikes          - 
%   plottwogrid             - 
%   PlotXcorrC              - Plots cross-correlations for a Clusters set
%   png2pgm                 - 
%   PrintCells              - PrintCells(C)
%   PrintErrors             - 
%   PrintMsg                - s = PrintMsg(logfid,format, args)   calls sprintf(format,args)
%   ProcessGridFullV        - DATA = ProcessGridFullV(name, varargin) make FullV files from .ns5 files. 
%   ProfileSummary          - ProfileSummary(P, varargin)  summarizes results from P = profile('info')
%   psfgui                  - 
%   psychdatelist           - 
%   PsychMon                - X = PsychMon(filename)
%   PsychSum                - ps = PsychSum(name, varargin)
%   PulseTrigLFP            - 
%   r2d                     - 
%   RadialSum               - sor = RadialSum(kernel, fixr)
%   Rayleigh                - 
%   RC2mat                  - x = RC2mat(rc, varargin)
%   ReadAutoLog             - 
%   ReadCellList            - read a list of .mat files for a list GUI
%   ReadConfig              - Read a file that sets the configuation of GUI parameters
%   ReadDXseq               - get dx: lines from serial output and match up with id numbers
%   ReadErrors              - print out errors for Data Folders
%   ReadExptDir             - [Expts, Idx] = ReadExptDir(name, varargin)
%   ReadExptFiles           - 
%   ReadHeader              - 
%   ReadManualExpt          - 
%   ReadOnlineTxt           - read online data with just spike times in text format
%   ReadPen                 - pen = ReadPen(file, varargin) Reads/Plots a penetration log file
%   ReadPGM                 - [im, details] = ReadPGM(name... read pgm image and comments
%   readpsych               - 
%   ReadPsychFile           - File format
%   readpsychsum            - 
%   ReadSerial              - 
%   ReadSerialFile          - [result, Expt, txt] = ReadSerialFile(name, 'readexpt')
%   ReadSerialLine          - 
%   ReadSerialTrials        - [Trials, txt] = ReadSerialTrials(name) Read serial file fo AplaySpkfile
%   ReadSpikeFile           - Spikes = ReadSpikeFile(spkfile, varargin)
%   ReadSpkFile             - 
%   ReadTemplateClusters    - ReadTemplateClustes(a)
%   ReadUfl                 - res =ReadUfl(name)
%   RebuildImages           - RebuildImages(Expt) Builds images for an ORBW Expt
%   RecalcSubspace          - [sdf, count] = RecalcSubspace(E, tx, xvals,trange) %recalculate an sdf.
%   Recombine               - Recombine(name,varargin) calls combine to reebuild name
%   ReCountExptSpikes       - 
%   RectAdd                 - M = RectAdd(a,b,...) creates a matrix M from adding
%   RemakeClusters          - fixed - RemakeClusters(D, varargin) Calls AllVPcs for Clusters that need
%   RemapSquareMatrix       - remapSquareMatrix(dm, map)  re-order a square matix using order defined in map
%   RemoveFullVMains        - 
%   ReplaySpikes            - 
%   resampler               - a is a vector of .mat psychfiles for ORBWCRC expts.
%   RestoreFile             - RestoreFile(name, ...    Restore Backups made with BackupFile
%   ReWriteClusters         - 
%   RewriteProbeSpikes      - 
%   RlsRCc                  - allframes = RlsRCa(Expt, seedlist, idlist, varargin)
%   rmfields                - rmfields(S, a,b,c....)    removes a list of fields from S
%   rotate2d                - 
%   rotatesac               - 
%   RunAllGridFiles         - RunAllGridFiles(file, .....
%   RunAllVPcs              - res = RunAllVPcs(dir, .....
%   runbinoc                - 
%   RunExptCmds             - res =  RunExptCmds(name,expts, varargin)
%   runlist                 - Basic list -running GUI.
%   RunPlotClusters         - 
%   SaveArrayConfig         - SaveArrayConfig(DATA) save array configuration
%   SaveCallback            - 
%   SaveComments            - 
%   SaveConfig              - SaveConfig(DATA, file, savefields, varargin)    
%   SaveFields              - SaveFields(name, X) save the fields of X in file name
%   SaveLayout              - DATA = ApplyLayout(DATA, varargin)
%   ScaleLFP                - 
%   ScaleLFPRC              - 
%   ScaleWindow             - resize a window and its contents
%   scanlines               - txt = scanlines(name)  retunrs a cell array of strings corresponding
%   sd2cv                   - converts sd in degrees to vector length (1-circular variance)
%   sdflatency              - [latency, details] = sdflatency(sdf, presamples, varargin) estmate latency
%   SelectFits              - gid = SelectFits(fits,varargin)
%   SelectProbe             - selected = SelectProbe(probelist, selected)
%   sem                     - r = sem(x) std(x)./sqrt(length(x))
%   Serial2Netexp           - [Trials, txt] = Serial2Netexp(name) Read serial file, and create
%   Serial2Stim             - ReadSeeds  gets seed/id numbers from online text record to fix missing
%   SerialTest              - 
%   SerialToAddTxt          - 
%   SetCheck                - SetCheck(tag, value)
%   SetClusters             - modifies ClusterTimesFiles on disk
%   SetData                 - SetData(DATA) shortcut for set(DATA.toplevel,'UserData',DATA);
%   SetDiag                 - 
%   SetExptRC               - Expt = SetExptRC(Expt, varargin) Make sure RC fields are set in Exot
%   SetExptSpikes           - [DATA, ispk, dprime, details] = SetExptSpikes(DATA, expid, show, varargin)
%   SetFigData              - 
%   SetFigPos               - SetFigPos(X, tag)
%   SetFigureName           - 
%   SetMenu                 - 
%   SetMenuCheck            - Turns on/off checked proerty in menu itesm
%   SetPdir                 - 
%   SetPlottingDefaults     - 
%   SetProbesToUse          - 
%   SetPropertyRecursive    - SetPropertyRecursive(X, prop, ...
%   SetSpkCodes             - 
%   SetSubplots             - SetSubplots(nrow, ncol, plots, flag)
%   SetTrialOffsets         - SetTrialOffsets(Expts, varargin) ensure Trials.Trial increases
%   SetTrialVals            - 
%   SetUIFont               - SetUIFont(F, varargin) sets the font for UI elements belonging to a figure
%   SetUIRecursive          - 
%   ShowExptErrs            - ShowExptErrs(Expt, varargin) Prints out errors in struct Expt 
%   ShowSlices              - 
%   ShuffleTimes            - st = ShuffleTimes(t, Expt,varargin)
%   simplegui               - Bareboes gui window with
%   sine                    - sine(params,x) make generate sinewave over x
%   sixteen2eight           - 
%   smhist                  - [smoothed, x] = smhist(dat, ...) make smoothed histogram
%   smooth                  - smooth(x, w, ...) smooth data with boxcar, width w
%   smooth2                 - 
%   smrcat                  - smrcat(outf, file1, file2, ...) Combines contents of .mat files made by Spike2
%   Spike2                  - Matlab files for spike 2
%   SpikeCycle              - [spikes, nt]  = SpikeCycle(Trials, tmin, tmax, period)
%   SpikeHist               - [w, dp] = SpikeHist(list, ...)
%   SpikeISI                - 
%   SpikeShape              - dp = SpikeShape(SPK, varargin)
%   Spk2Double              - Spk2Doulbe(S)Converts ints to double in Spks structs
%   SpkDefs                 - 
%   Spks2Sdf                - 
%   split                   - cellstr = split(s,delimiter)  splits a string into a cell array
%   SplitMatlab             - SplitMatlab(name, @class) Splits a .m file into separate files for each function
%   splitname               - [dir, name, fitname] = Splitname(path)
%   splitpath               - [name, dir] = splitpath(path, varargin) replicate fileparts
%   SpoolFullV              - AllVPcs(V, ...)  takes an MxN matrix (electrode by voltage) of continuous
%   SpoolSpikes             - GUI for playing back spikes from existing expts.
%   sprintcell              - s = sprintcell(X, varargin) build a string from a cell array
%   SpTrigLFP               - SpTrigLFP(Trials, duration, samplerate, w)
%   SqrtGaborFitSSD2        - 
%   squishDistanceMatrix    - squishDistanceMatrix rearranges the input distance matrix to maximize the 
%   startup                 - 
%   StimulusName            - Convert stimulus name string to num or vice versa
%   StoreData               - 
%   strstr                  - [found, ids] = strstr(str, pattern) 
%   StrTable                - print a cell array of strings as a table
%   subsample               - subsmp = subsample(x, ratio)
%   SummarizeCells          - 
%   Swatch2Clusters         - Swatch2Clusters(name, varargin) Convert old files to PlotClusters
%   SxCx                    - caluclate modulation in Autocorrelatin functions to classify Simple/Compex
%   tab2im                  - 
%   Table2Expt              - 
%   TestEllipse             - 
%   TestEllipseDraw         - 
%   TestEvents              - 
%   TestFullV               - 
%   testguide               - This is the machine-generated representation of a Handle Graphics object
%   testscan                - Experiment with methods for reading psych files and other text data
%   TestSpkFile             - [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%   TestSpkFileA            - [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%   TestSpkFileB            - [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%   TestStimOn              - TestStimOn
%   TestText                - 
%   TestTrigsdf             - 
%   tilefigs                - <cpp> tile figure windows usage: tilefigs ([nrows ncols],border_in pixels)
%   timediff                - timediff(t, varargin) calculates differences between time values, in seconds (default)
%   TimeMark                - tt = TimeMark(tt, str, show)  keep track of exection time
%   TimeRange               - 
%   TimesInTrial            - TimesInTrial(t, Expt, preperiod, postperiod)
%   TrackBlanks             - 
%   TrackFigPos             - called when fiugres are closed to store current position
%   TrackLFPMains           - TrackLFPMains  adjusts an  LFP file for flucutating mains noise.
%   TreeFind                - [names, sizes, dates, pathnames] = TreeFind(path, varargin)
%   TreeSummary             - [res, details] = TreeSummary(path, varargin) Summarize disk usage
%   TrigLFP                 - TrigLFP(Trials, times, samplerate, nch)
%   trigsdf                 - [sdf, nsac, nspikes, spikes, times] = trigsdf(Trials, width, times, ...)
%   trigsdfa                - [sdf, nsac, nspikes, spikes, times] = trigsdf(Trials, width, times, ...)
%   trigsdfb                - [sdf, nsac, nspikes, spikes, times] = trigsdf(Trials, width, times, ...)
%   trigsdfc                - [sdf, nsac, nspikes, spikes, times] = trigsdf(Trials, width, times, ...)
%   TryFullV                - CheckFullV(FullV, varargin)
%   twopoolcounts           - [acounts, bcounts, details]=twopoolcounts(r,rb,ntrials,poolsize,varargin)
%   TwoProbes               - [aid, bid] = TwoProbes(sa,sb, varargin)
%   untitled2               - 
%   UpdateBabProbeList      - FInds probes markded as bad in Clusters, and Adds these to the ArrayConfig
%   UpdatePsychFile         - UpdatePsychFile(name, copies local pysch file to network
%   UpdateSpike2Files       - 
%   urndwalks               - 
%   VallTimes2id            - id = VallTimes2id(Vall, t) 
%   verg                    - binoc
%   VergPsych               - standalone version of Psych plotter from verg.
%   vergverion              - 
%   vergversion             - 
%   WeightedSum             - [m, n] = WeightedSum(x, n, varargin)
%   WinFront                - WinFront(tag)
%   winstartup              - 
%   WorkerString            - Return string with worker number 
%   xcorrtimes              - [xc, details] = xcorrtimes(a, b, varargin)
%   xy2op                   - xy2op(op, angle,..)
%   xy2pol                  - convert Cluster.xy to radius/angle
%   XYClassify              - 
%   xyrotate                - XY = xyrotate(x,y,angle)
%   xysmooth                - [sx,xy] = xsmooth(x, y, w, ...) smooth data with boxcar, width w
