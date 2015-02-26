% @ALLV
%
% Files
%   ._AllVPcs               - 
%   ._OriginalAllVPcs       - 
%   AddAxisContextMenu      - 
%   AddCellLabel            - 
%   AddCellMean             - 
%   AddCellMenu             - 
%   AddComment              - 
%   AddErr                  - 
%   AddLineContextMenu      - 
%   AddMarkToPCs            - 
%   AddParameterMenu        - 
%   AddSelectorContextMenu  - 
%   AllV                    - 
%   AllVPcs                 - AllVPcs(V, ...)  takes an MxN matrix (electrode by voltage) of continuous
%   ApplyLayout             - 
%   AssignCluster           - 
%   AutoCut                 - 
%   AutoCutAll              - 
%   AutoCutOne              - 
%   BadCluster              - 
%   BestAngle               - Find best angle to cut in a 2D space. With very skewed distributison
%   BestAngleGM             - Find best angle to cut in a 2D space, using 1-D Gaussian Mixtures
%   BestGMAngle             - 
%   BestSpace               - 
%   BimodalCoeff            - 
%   BoundaryFromCluster     - 
%   BuildAllV               - 
%   ButtonPressed           - 
%   ButtonReleased          - 
%   CalcClusterDistance     - 
%   CalcDistanceMatrices    - overfit with GMs, and calculate distance Matrix
%   CalcDprime              - 
%   CalcICA                 - 
%   CalcISI                 - [isis, trials, spkids] = AllV.CalcISI(Trials, varargin)
%   CalcPCs                 - 
%   CalcRadius              - 
%   CalcScores              - this need to be modified to work like AllV.TemplatePlot now, with
%   CalcTemplateScores      - 
%   CalcTemplatesFromMean   - 
%   CalculateTemplateDips   - 
%   CalcXcorr               - 
%   Cell2Cluster            - 
%   ChangeCell              - 
%   ChangeProbe             - 
%   CheckClusterFields      - 
%   CheckClusterMarks       - 
%   CheckClusters           - 
%   CheckClusterValues      - 
%   CheckForMean            - 
%   CheckScoreScaling       - Template scores are all normalized so that spaces are roughly isotropic
%   CheckSign               - 
%   CheckTemplates          - 
%   ClassifyAll             - 
%   ClassifyAndFit          - 
%   ClassifyFit             - 
%   ClassifySpikes          - 
%   CloseLog                - 
%   ClusterFile             - 
%   ClusterFromBoundary     - 
%   ClusterFromPoints       - 
%   ClusterInfo             - 
%   ClusterIsSet            - check that cluster cl is defined for a give clusterstruct
%   CompareMean             - 
%   CondenseCluster         - 
%   CondenseClusters        - 
%   CutAndPlot              - 
%   CutAndSave              - 
%   DistanceToEllipse       - When called from a new graph, ther is nothing to check, but E may have
%   DrawEllipse             - 
%   DrawLine                - 
%   EllipseDip              - [SSD, details ] = AllV.MinimiseEllipseb(params, DATA, state)
%   ExcludeTrials           - 
%   exitallv                - 
%   ExptFigMenu             - 
%   FindNearestCluster      - 
%   FinishSpikePlot         - 
%   Fit2Gauss               - [dp, fits, details] = AllV.Fit2Gauss(C, r, DATA, varargin)
%   FitGaussMeans           - 
%   FitWindow               - 
%   FullVKeyPressed         - 
%   GetAllV                 - 
%   GetClusterDef           - return relevant member of teh clsuter struc    
%   GetClusterPlot          - 
%   GetClusterPlots         - 
%   GetClusterXYData        - 
%   GetCSD                  - 
%   GetDataStruct           - 
%   GetGuiState             - 
%   GetProbeFromName        - 
%   GetSubCluster           - return struct with details for a given cluster
%   GetValues               - 
%   gmdistance              - 
%   gmdprime                - calcualte drpime between two Gaussians in gmdistribution fit        
%   GMfit                   - 
%   GMMButton               - 
%   GMMmenu                 - 
%   GoodCluster             - 
%   HistButtonDragged       - 
%   HistButtonPressed       - 
%   HistButtonReleased      - 
%   HistKeyPressed          - 
%   HistMenu                - 
%   HitImage                - 
%   HitISI                  - 
%   HitTrial                - 
%   HitXYPlot               - 
%   IDStr                   - 
%   InGraph                 - 
%   isacell                 - [true, cellid] = AllV.isacell(DATA, ei, p)
%   IsTemplateCut           - 
%   IterateFit              - 
%   IterateTemplateFit      - 
%   JamesAutoCut            - 
%   KeyPressed              - 
%   LineLength              - 
%   LoadCellFile            - 
%   LoadClusterDetails      - 
%   LoadComments            - 
%   LoadDataClusters        - 
%   LoadExpt                - [Expt, matfile] = AllV.LoadExpt(DATA, ei)
%   LoadExptA               - 
%   LoadTrigTimes           - 
%   MakeJamesSpikes         - 
%   MarkAxes                - 
%   MaxDprime               - 
%   memsize                 - 
%   MinimiseEllipse         - [SSD, details ] = AllV.MinimiseEllipse(params, DATA, state)
%   MinimiseEllipseb        - [SSD, details ] = AllV.MinimiseEllipseb(params, DATA, state)
%   MinimiseEllipseDip      - 
%   MinimiseLine            - 
%   MiscMenu                - 
%   mygetappdata            - 
%   myisappdata             - 
%   mysetappdata            - 
%   NeedClusterData         - 
%   NeedMore                - Cluster parameters can look bad if ALL of the triggered events are from a
%   NeedTemplateForCluster  - 
%   NextPCs                 - 
%   oldDrawEllipse          - 
%   oldDrawLine             - 
%   oldFindDip              - 
%   OldPlotFullV            - 
%   OldSetMenuCheck         - 
%   OptimizeBoundary        - 
%   OptimizeClusterBoundary - 
%   OptimizeEllipse         - 
%   OptimizeLine            - 
%   OptimizeVarE            - 
%   OptionMenu              - 
%   OriginalAllVPcs         - AllVPcs(V, ...)  takes an MxN matrix (electrode by voltage) of continuous
%   PCButtonDragged         - 
%   PCCluster               - 
%   PCKeyPressed            - 
%   PCLabels                - 
%   PlaceUi                 - 
%   Plot2GaussFit           - 
%   PlotAllMeans            - 
%   PlotAllProbes           - 
%   PlotAllXCorr            - 
%   PlotCellList            - 
%   PlotCluster             - --- help for PC/PlotCluster ---
%   PlotClusters            - --- help for PC/PlotClusters ---
%   PlotClusterXY           - 
%   PlotExptCounts          - 
%   PlotFeatures            - 
%   PlotFullV               - 
%   PlotGauss2Fit           - 
%   PlotGridSpikes          - 
%   PlotHistogram           - 
%   PlotISI                 - 
%   PlotMeanSpike           - 
%   PlotMeanSpikes          - 
%   PlotMenu                - 
%   PlotOneTrial            - 
%   PlotOneXY               - 
%   PlotPCs                 - 
%   PlotProbeSpikes         - 
%   PlotQuickSpikes         - 
%   PlotResult              - 
%   PlotSpikes              - 
%   PlotTemplateScores      - 
%   PlotTrialSpikes         - 
%   PlotTriggerHist         - 
%   PlotVals                - 
%   PlotVarE                - 
%   PlotXcorr               - 
%   PlotXY                  - 
%   PlotXYSequence          - 
%   PreferencesPopup        - 
%   PrevCluster             - 
%   ProbeMenu               - 
%   ProbeNumber             - AllV.ProbeNumber returns Real Probe number
%   ProbeSelector           - 
%   ProjectND               - 
%   QuantifyQuickClusters   - 
%   QuickAutoCut            - 
%   QuickSpks               - 
%   ReadFromLog             - 
%   ReadSpikeFiles          - 
%   ReClassify              - 
%   ReplotPCs               - 
%   ReplotXcorrs            - 
%   ResetDataForNewProbe    - 
%   RestrictTimeRange       - 
%   RetriggerDialog         - 
%   RotateCluster           - just in case its called before cluster is set
%   RotateLine              - 
%   Rprime                  - 
%   RunGMMFit               - 
%   SaveCellList            - 
%   SaveClusters            - 
%   SaveMeanSpikeOnly       - 
%   SaveSpikes              - 
%   ScrollSpikes            - 
%   ScrollV                 - 
%   ScrollWheel             - 
%   SelectProbe             - 
%   SelectTrial             - 
%   SetADC                  - 
%   SetCellCompare          - 
%   SetCellEntry            - 
%   SetCellFromLine         - 
%   SetClusterIcon          - 
%   SetEllipseDrawing       - 
%   SetFigure               - [F, isnew] = AllV.SetFigure(tag, DATA, ..)
%   SetFullVNames           - 
%   SetGUI                  - 
%   SetMenuChecks           - Set checked on/off for menu items in hm whose tags
%   SetOption               - 
%   SetPCs                  - 
%   SetPlot                 - 
%   SetTemplateData         - 
%   SetTrialList            - 
%   SetupVall               - 
%   SetVsamples             - DATA.vsmps = [20 6 15 11 30 20];
%   ShapeCorr               - 
%   ShowADCPos              - 
%   ShowFullV               - 
%   ShowTaggedProbes        - 
%   SmallCluster            - remove fields from C that use memory
%   SpikeButtonPressed      - 
%   SpikeDraw               - 
%   SpkFileName             - 
%   SpoolAllSpikes          - 
%   SpoolSpikes             - 
%   StripClusters           - 
%   SummaryHit              - 
%   SyncSpikes              - 
%   TemplateGMFits          - 
%   TemplateLabels          - 
%   TemplatePlot            - 
%   TemplateSpace           - 
%   Text2Val                - 
%   TimeMark                - 
%   TrialMarkChar           - 
%   TriggerV                - 
%   UseAllEvents            - 
%   UseProbeList            - 
%   WhichPlotType           - 
%   XYplot                  - 
%   xyr2pos                 - 
