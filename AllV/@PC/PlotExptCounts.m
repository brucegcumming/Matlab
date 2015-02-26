 function Expt = PlotExptCounts(DATA, e, p, cl, varargin)     Expt = [];     if strncmpi(DATA.plotexpttype,'none',4)           if isfield(DATA,'Expt')             Expt = DATA.Expt;         end         return;     end     if nargin < 4         cl = DATA.currentcluster;     end     if nargin < 3        p = DATA.currentpoint(2);     end     if nargin < 2        e = DATA.currentpoint(1);     end     tag = DATA.tag.expt; %set to empty to not set figure       plotargs = {};       j = 1;       while j <= length(varargin)           if strcmpi(varargin{j},'hold')               plotargs = {plotargs{:} 'hold'};           elseif strcmpi(varargin{j},'tag')               j = j+1;               tag = varargin{j};           end           j = j+1;       end                     if ~isempty(tag)           PC.SetFigure(DATA,tag);       end       eid = DATA.exptid(e);     [Clusters, DATA] = PC.CheckClusterLoaded(DATA, e);     Expts = getappdata(DATA.toplevel,'Expts');     clnum = cl;     Expt = PC.CountExptSpikes(DATA, Expts{eid},Clusters{e}{p},clnum);       if strncmp(Expt.Header.expname,'image.orXob',11) % Plotexpt does this right now, for psych and non-psych types           plotargs = {plotargs{:}};       end     if strncmpi(DATA.plotexpttype,'trialcounts',10)         PlotExpt(Expt,'seqt','rcnmin',10,plotargs{:});     elseif strncmpi(DATA.plotexpttype,'none',4)     else         PlotExpt(Expt,'shown','rcnmin',10,'fbox',plotargs{:});     end