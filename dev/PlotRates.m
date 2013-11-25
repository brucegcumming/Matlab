function [result, Expt] = PlotRates(Expt,type, varargin)

% result = PlotRates(Expt,type, ...)
% divides trials from Expt into groups according to value of the
% field specified by 'type', calcuales mean and SD of these, then
% plots these mean rates as a function of parameter 'type'
%
% returns result, a structure with four vectors, whose length is
% the number of different stimuli. 
% result.x        the value of the stimulus parameter.
% result.means    the mean spike countg
% result.n        number of trials in the average.
% result.sd       S.D. of the rate.
% result.flags    always zero (!!)
%
% PlotRates(Expt,type, 'Type2', type2) plots several lines, one for
% each unique value of the field 'type2'
%
% PlotRates(Expt,type, 'LogX') used log X axes.
%
% PlotRates(Expt,type, 'RandPhase') separates out any trials with
% random phase and plots these as an "extra" interleaved condition
%
% PlotRates(Expt,type, 'Extra', 'expr', 'label')  separates out any
% trials that satisfy    eval('[Expt.Trials.' expr);
% and plots these separately, with the label 'label'
% e.g. PlotRates(Expt, 'dx', 'Extra', '.ce == 0', 'UnCorr')
% will plot uncorrelated trials separately, with the label UnCorr
% in the legend.
%
% PlotRates(E%pt,type, 'sdf') makes sdfs
% PlotRates(Expt,type, 'sdf','showsum') adds sdf avg across all trials
%                      'sdfall')  does the same
%
%PlotRates(Expt, type, 'sval','code',x)
% selects a subset of trials where Expt.Trials.(code) == x
% PlotRates(Expt,type, ) 


gettimes =1;
colors = mycolors;
fontsize = 12;
logx = 0;
addn = 0;
addl = 0;
plotlfp = 0;
legendpos = 0;
times = [];
mksdf = 0;
sdfw = 100;
showrates = 0;
fillsymbols = 1;
fillall = 0;
laxis = 0;
raxis = 0;
getmod = 0;
showpcolor = 0;
forcezero = 1;
showtitle = 1;
showplot = 1;
showsum = 0;
verbose = 0;
dosqrt = 0;
autonmin = 0;
triglfp =0;
psychtag = 'Psych';
symbols = 'ooooooooooooooooooooooooooooooooo';
lfprange = [30 60 200 10 4];
nlfp = 1;
linestyle = {'-',':','-.','--','-', ':', '-.', '--' };
nline = 1;
nolines = 0;
acov = 0;
plotsimple = 0;
showmu = -1;
sethold = 0;
showsdf = []; 
subsdf = '';
preperiod = 2000;
postperiod = 2000;
includebad = 0; %set to 1 to include invalid saccade trials for tuning curves
result = [];
colorids = [];


if isstruct(type)
    result = type;
    type = result.type{1};
end

if(~exist('Expt') | ~isfield(Expt.Trials,type))
    result = [];
  return;
end

extraexp{1} = 'st] == 0';
extra.label{1} = 'Blank';
nextra = 1;
splitextra(1) = 0;
showcounts = 0;
periodsdf = 10;
periodic = 0;
sdfargs = {};
psfargs = {};
duration = [];
nflp = 1;
latency = 500;
forcexvs = 0;
xvs = [];

nvar = nargin - 2;
nmin = 1;
j = 1;
 while j  <= nvar
    str = varargin{j};
    if strncmpi(str,'Type2',5)
        j = j+1;
        if isempty(strmatch(varargin{j},{'e0' type}))
            btype = varargin{j};
            if isfield(Expt.Trials,btype)
                bvals = eval(['sort(unique([Expt.Trials.' btype ']));']);
            else
                clear btype;
            end
        else
            clear btype;
        end
    elseif strncmpi(str,'Type3',5)
        j = j + 1;
        ctype = varargin{j};
        if isfield(Expt.Trials,ctype)
            cvals = eval(['sort(unique([Expt.Trials.' ctype ']));']);
        else
            clear ctype;
        end
    elseif strncmpi(str,'Type',4)
        j = j + 1;
        type = varargin{j};
    elseif strncmpi(str,'showmu',5)
        showmu = 0;
    elseif strncmpi(str,'block',5)
      j = j+1;
      block = varargin{j};
      varargin = {varargin{[1:j-2 j+1:end]}};
      j = j -2;
      nvar = nvar -2;
      Trials = [];
      for k = 1:length(block)
          if block(k) < length(Expt.Header.BlockStart)
              tid = Expt.Header.BlockStart(block(k)):Expt.Header.BlockStart(block(k)+1)-1;
              id = find(ismember([Expt.Trials.Trial],tid));
              if length(id)
              Trials = [Trials Expt.Trials(id)];
              else
                  fprintf('No Trials for Block %d\n',block(k));
              end
          end
      end
      if length(Trials) > 2
          Expt.Trials = Trials;
      end
    elseif strncmpi(str,'acov',4) %calculate autocorrelation
        acov = 1;
    elseif strncmpi(str,'cmpduration',6)
      PlotRates(Expt,type,varargin{1:j-1});
      hold on;
      j = j+1;
      duration = varargin{j};
    elseif strncmpi(str,'colorids',6)
      j = j+1;
      colorids = varargin{j};
    elseif strncmpi(str,'colors',6)
      j = j+1;
      colors = varargin{j};
  elseif strncmpi(str,'duration',4)
      j = j+1;
      duration = varargin{j};
  elseif strncmpi(str,'fillall',5)
      fillall = 1;
    elseif strncmpi(str,'forcecolor',8)
      j = j+1;
      forcecolors = 1;
      for k = 1:50
      colors{k} = varargin{j};
      end
      symbols = 'os+xos+xos+xos+xos+xos+xos+xos+x';
  elseif strncmpi(str,'Hold',4)
      hold on;
      sethold = 1;
  elseif strncmpi(str,'latency',4)
      j = j+1;
      latency = varargin{j};
  elseif strncmpi(str,'line',4)
      j = j+1;
      nline = varargin{j};
  elseif strncmpi(str,'Mod',3)
      getmod = 1;
  elseif strncmpi(str,'periodic',5)
      periodic = 1;
  elseif strncmpi(str,'sdperiodic',5)
      j = j+1;
      periodsdf = varargin{j};
  elseif strmatch(str,'LogX')
    nextra = nextra+1;
    extraexp{nextra} = [type '] == 0'];
    extra.label{nextra} = 'Zero';
    splitextra(nextra) = 0;
    logx = 1;
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
elseif strncmpi(str,'Noplot',5)
      showplot = 0;
  elseif strncmpi(str,'NoFill',5)
      fillsymbols = 0;
  elseif strncmpi(str,'Noline',5)
      linestyle = {'none','none'};
  elseif strncmpi(str,'NoZero',5)
      forcezero = 0;
  elseif strmatch(str,'Nlin')
      addn = varargin{j+1};
      j = j+1;
  elseif strncmpi(str,'Nmin',4)
      if ischar(varargin{j+1}) & strncmpi(varargin{j+1},'auto',4)
          autonmin = 1;
      else
          nmin = varargin{j+1};
      end
      j = j+1;
  elseif strmatch(str,'legendpos')
      legendpos = varargin{j+1};
      j = j+1;
  elseif strcmpi(str,'Extra')
    nextra = nextra+1;
    extraexp{nextra} = varargin{j+1};
    extra.label{nextra} = varargin{j+2};
    splitextra(nextra) = 0;
    j = j+2;
elseif strncmpi(str,'pcolor',5)
    showpcolor = 1;
    fillall = 1;
elseif strncmpi(str,'plotall',5)
    showplot = 2;
elseif strncmpi(str,'preperiod',8)
    j = j+1;
    preperiod = varargin{j};
elseif strncmpi(str,'postperiod',8)
    j = j+1;
    postperiod = varargin{j};
elseif strncmpi(str,'figPsych',4)
      psychtag = varargin{j+1};
      j = j+1;
elseif strncmpi(str,'Psych',4)
    if ~isfield(Expt.Trials,'RespDir')
        result = [];
          return;
    end
    if ~isfield(Expt.Trials,'Id')
        for t = 1:length(Expt.Trials)
            Expt.Trials(t).id = Expt.Trials(t).Trial;
        end
    end
%RespDir -1 = upward saccade (vs positive), = FAR choice for disparity task  
%for cylinder task this make a mess meaning of RespDir changes at 90/-90
% call with 'psychflip'
      if strncmpi(str,'Psychflip',8)
          uidx = find([Expt.Trials.RespDir] > 0);
          idx = find([Expt.Trials.RespDir] < 0);
      else
          uidx = find([Expt.Trials.RespDir] == -1);
          idx = find([Expt.Trials.RespDir] == 1);
      end

      if j > 1
          args = {varargin{1:j-1}};
          if j < length(varargin)
              args = {args{:} varargin{j+1:end}};
              k = j+1;
              while k < length(varargin)
                  if strncmpi(varargin(k),'hold',4)
                      sethold = 1;
                  end
                  k = k+1;
              end
          end
      else
          args = {};
      end
      %make sure same duration is calculated for both choices
      if isempty(duration)
              duration = min([Expt.Trials([uidx idx]).End] - [Expt.Trials([uidx idx]).Start]);
              args = {args{:} 'duration' duration};
      end
      if strncmpi(str,'Psychnoplot',8)
          args = {args{:} {'noplot'}};
          showplot = 0;
      end
      xvs = sort(unique([Expt.Trials.(type)]));
      xvs = xvs(find(~isnan(xvs)));
      uExpt = Expt;
      uExpt.Trials = Expt.Trials(uidx);
      dExpt = Expt;
      dExpt.Trials = Expt.Trials(idx);
      if ~isempty(uidx)
          [ures, exp] = PlotRates(uExpt,type,'fillall','xval',xvs,args{:});
          for tr = 1:length(uidx)
              Expt.Trials(uidx(tr)).count = exp.Trials(tr).count;
          end
          ures.RespDir = mean([uExpt.Trials.RespDir]);
          if showplot
              hold on;
          end
      end
      if ~isempty(idx) & ~isempty(uidx)
          [dres, exp] = PlotRates(dExpt,type,'Nofill','fillall','xval',xvs,args{:},'line',2,'hold','nmin',1);
          errstr = [];
          for tr = 1:length(idx)
              Expt.Trials(idx(tr)).count = exp.Trials(tr).count;
          end
          npvals = min([size(dres.counts,1) size(ures.counts,1)]);
          dres.RespDir = mean([dExpt.Trials.RespDir]);
          if sum(size(dres.x) ~= size(ures.x)) > 0
          [ys, ya, yb] = intersect(dres.y(1,:),ures.y(1,:));
          [xs, xa, xb] = intersect(dres.x(:,1),ures.x(:,1));
          dres.x = dres.x(xa,ya);
          dres.y = dres.y(xa,ya);
          dres.n = dres.n(xa,ya);
          dres.linevals = dres.linevals(ya);
          dres.colors = dres.colors(ya);
          dres.means = dres.means(xa,ya);
          dres.counts = dres.counts(xa,ya);
          dres.ids = dres.ids(xa,ya);
          if max(ya) <= length(dres.handles)
          dres.handles = dres.handles(ya);
          end
          ures.x = ures.x(xb,yb);
          ures.y = ures.y(xb,yb);
          ures.n = ures.n(xb,yb);
          ures.means = ures.means(xb,yb);
          ures.ids = ures.ids(xb,yb);
          ures.counts = ures.counts(xb,yb);
          if max(yb) <= length(ures.handles)
          ures.handles = ures.handles(yb);
          end
          ures.linevals = ures.linevals(yb);
          ures.colors = ures.colors(yb);
          if isfield(dres,'sdfs')
              dres.sdfs = dres.sdfs(xa,ya);
              ures.sdfs = ures.sdfs(xb,yb);
          end
          if isfield(dres,'labels')
              dres.labels = dres.labels(xa,ya);
              ures.labels = ures.labels(xb,yb);
          end
          end
          nx = 0;
          plotflip = 0;
          for k = 1:npvals;
              for ij = 1:min([size(dres.counts,2) size(ures.counts,2)])
                  dres.psum(k,ij) = length(dres.counts{k,ij}) + length(ures.counts{k,ij});
                  dres.presp(k,ij) = length(dres.counts{k,ij});
                  if strcmp(Expt.Stimvals.e2,'Id') || strcmp(Expt.Stimvals.e2,'e0') || strcmp(type,'cvsign')
                      dres.psychval(k,ij) = dres.x(k,ij);
                  else
                      dres.psychval(k,ij) = dres.y(k,ij);
                  end
                  nx = nx+1;
              end
          end
          if nx == 0
              dres.psum = 0;
              dres.presp = 0;
              dres.psychval = [];
          end
          ures.psum = dres.psum;
          ures.presp = dres.presp;
          ures.psychval = dres.x;
          if showplot
              figid = gcf;
              GetFigure(psychtag);
              if sethold == 0
              hold off;
              end
              pres = ExptPsych(Expt,'forcecolor',colors{1});
              showplot = 2;
          end
          if strcmp(Expt.Stimvals.e2,'Id') || strcmp(Expt.Stimvals.e2,'e0')
              for k = 1:size(dres.psum,2)
                  ps = dres.presp(:,k)./dres.psum(:,k);
                  h(k) = plot(dres.psychval(:,k),ps,'o-','color',colors{k});
                  labels{k} = sprintf('%.2f',dres.linevals(k));
                  if showcounts
                      step = range(dres.psychval(:))/20;
                      for ti = 1:size(dres.psychval,1)
                          text(dres.psychval(ti,k)+step,ps(ti),sprintf('%d',dres.psum(ti,k)),'color',colors{k});
                      end
                  end
                  hold on;
              end
              legend(h,labels);
          elseif size(dres.psum,2) == 2
              
              if sum(size(dres.counts)) == sum(size(ures.counts))
              for k = 1:size(dres.counts,2)
                  for ij = 1:npvals
                      dres.psum(ij,k) = length(dres.counts{ij,k}) + length(ures.counts{ij,k});
                      dres.presp(ij,k) = length(dres.counts{ij,k});
                      dres.psychval(ij,k) = dres.x(ij,k);
                  end
              end
              else
                  for k = 1:size(dres.counts,2)
                      for ij = 1:npvals
                          dres.psum(ij,k) = length(dres.counts{ij,k});
                          dres.presp(ij,k) = length(dres.counts{ij,k});
                          dres.psychval(ij,k) = dres.x(ij,k);
                      end
                  end
              end
              len = size(dres.psum,1);
              %zo inverts the scale for Ori Bandwidtih, since large bandwidths = 0 signal
              if strmatch(btype,'ob')
                  zo = dres.psychval(len,1);
                  X = sd2cv(dres.psychval) .* sign(dres.x-mean(dres.x(:)));
                  id = find([Expt.Trials.ob] > 120);
                  [a,details] = CalcConsistency(Expt.Trials(id));
                  if details.rptfraction > 0.1
                      errstr = sprintf('%d/%d Seeds NOT repeated',sum(details.srpts==1),length(details.seeds));
                  end
              else
                  if strmatch(type,{'ob'})
                      X = sd2cv(dres.psychval) .* sign(dres.y-mean(dres.y(:)));
                  elseif strmatch(type,{'Dc'})
                      X = (dres.psychval) .* sign(dres.y-mean(dres.y(:)));
                  else
                      X = (dres.psychval) .* sign(dres.x-mean(dres.x(:)));
                  end
                  zo = 0;
              end
              for k = 1:prod(size((X)))
                  psf(k).n = dres.psum(k);
              psf(k).resp = dres.presp(k);
              psf(k).x = X(k);
              end
              if zo
                  psf(len).n = dres.psum(len,1) + dres.psum(len,2);
                  psf(len).resp = dres.presp(len,1) + dres.presp(len,2);
                  psf(len).x = 0;
              end
              if length(psf) > 1
                  fit = fitpsf(psf);
                  if showplot == 1 %now plotted with ExptPsych
                      plot([psf.x],[psf.resp]./[psf.n],'o','color',colors{k});
                      h = fitpsf(fit.data,'showfit',fit,psfargs{:},'color',colors{1});
                      fit.fith = h.fith;
                      set(gca,'Ylim',[0 1]);
                      title(sprintf('Bias %.1f, Slope %.3f',fit.fit(1),abs(fit.fit(2))));
                  end
              dres.psych = fit;
              ures.psych = fit;
              else
                  dres.psych = [];
                  ures.psych = [];
              end
          elseif strcmp(Expt.Stimvals.e2,'ob') && ~strcmp(type,'cvsign')
              id = find([Expt.Trials.ob] > 120);
              if isfield(Expt.Trials,'se')
              [a,details] = CalcConsistency(Expt.Trials(id));
              if details.rptfraction > 0.1
                  errstr = sprintf('%d/%d Seeds NOT repeated',sum(details.srpts==1),length(details.seeds));
              end
              end
              if sum(size(dres.psychval) == size(dres.x)) == 2
              X = sd2cv(dres.psychval) .* sign(dres.x-mean(dres.x(:)));
              for j = 1:prod(size(X))
                  pp(j).x = X(j);
                  pp(j).n = dres.psum(j);
                  pp(j).resp = dres.presp(j);
              end
              fit = fitpsf(pp);
              dres.psych = fit;
              if showplot
              plot([pp.x],[pp.resp]./[pp.n],'o');
              if isfield(fit,'data')
              h = fitpsf(fit.data,'showfit',fit,psfargs{:},'color',colors{1});
              end
                      title(sprintf('%sBias %.1f, Slope %.3f',errstr,fit.fit(1),abs(fit.fit(2))));
              end
              else
                  dres.psych = [];
              end
              ures.psych = [];
%              plot(sd2cv(dres.psychval) .*sgn,dres.presp./dres.psum,'o');
          else
                  for j = 1:size(dres.x,1)
                      pp(j).x = dres.x(j);
                      pp(j).n = dres.psum(j);
                      pp(j).resp = dres.presp(j);
                  end
                  fit = fitpsf(pp);
                  dres.psych = fit;
                  ures.psych = fit;
              if showplot == 1
                  plot(dres.psychval,dres.presp./dres.psum,'o');
                  if isfield(fit,'data')
                      h = fitpsf(fit.data,'showfit',fit,psfargs{:},'color',colors{1});
                  end
              end
          end
          if showplot
              figure(figid);
          end
          clear result;
          result(1) = ures;
          result(2) = dres;
      elseif ~isempty(idx)
          [dres, exp] = PlotRates(dExpt,type,'Nofill','fillall','xval',xvs,args{:});      
      end
      
      if ~exist('result','var')
                result = [];
      end
      j = j+1;
      return;
  elseif strncmpi(str,'Sdfw',4)
      j = j+1;
      sdfw = varargin{j};
  elseif strncmpi(str,'Sacsdf',6)
      mksdf = 2;
  elseif strncmpi(str,'Sdf',3)
      mksdf = 1;
      if strncmpi(str,'Sdfall',6)
          showsum = 1;
      elseif strncmpi(str,'Sdfsub',6)
          j = j+1;
          subsdf = varargin{j};
      end
  elseif strncmpi(str,'sxcx',4) %calculate autocorrelation
        acov = 1;
        plotsimple = 1;
  elseif strncmpi(str,'notitle',5)
      showtitle = 0;
  elseif strncmpi(str,'showsum',7)
      showsum = 1;
  elseif strncmpi(str,'sqrt',4)
      dosqrt = 1;
  elseif strncmpi(str,'lfp',3)
      plotlfp = 1;
      len = Expt.Header.lfplen;
      if isfield(Expt.Header,'LFPsamplerate')
          ftfrq = (0:len-1)/(len * Expt.Header.LFPsamplerate .* 10);
      else
          ftfrq = (0:len-1)/(len/(10000 * Expt.Header.CRsamplerate));
      end
      result.lfpfrq = ftfrq;
      subplot(2,1,2);
      hold off;
      if strncmpi(str,'lfptrig',5)
          triglfp = 1;
      elseif strncmpi(str,'lfpt',4)
          plotlfp = 2;
      end
  elseif strncmpi(str,'Trialids',8)
      j = j+1;
      tid = find(ismember([Expt.Trials.id],varargin{j}));
      Expt.Trials = Expt.Trials(tid);
  elseif strncmpi(str,'Trials',3)
      j = j+1;
      Expt.Trials = Expt.Trials(varargin{j});
  elseif strncmpi(str,'Rate',3)
      showrates = 1; 
  elseif strmatch(str,'AddMonoc') & isfield(Expt.Trials,'me')
      if exist('btype','var') & strncmp(btype,'co',2)
          nb = length(bvals);
          covals = bvals;
      elseif isfield(Expt.Trials,'co')
          idx = find([Expt.Trials.me] ~= 0);
          covals = unique([Expt.Trials(idx).co]);
          nb = length(covals);
      else
          nb = 1;
      end
% special care is required when different contrast combinations are used.
% Need to plot each monocular/contrast combination.
%
      if isfield(Expt.Trials,'c0')
          clvals = unique([Expt.Trials(idx).c0]);
          crvals = unique([Expt.Trials(idx).c1]);
          nb = max([length(crvals) length(clvals)]);
      end
      if isfield(Expt.Trials,'me')
          for k = 1:nb
              nextra = nextra+1;
              if nb > 1 & isfield(Expt.Trials,'c0')
                  if k <= length(clvals)
                      lbl = sprintf(' co %.2f',clvals(k));
                      extraexp{nextra} = sprintf('me] == 1 & abs([Expt.Trials.c0] - %.2f) < 0.1',clvals(k));
                  end
                  if k <= length(crvals)
                      lbr = sprintf(' co %.2f',crvals(k));
                      extraexp{nextra+1} = sprintf('me] == -1 & abs([Expt.Trials.c1] - %.2f) < 0.1',crvals(k));
                  end
              else
              extraexp{nextra} = 'me] == 1';
              extraexp{nextra+1} = 'me] == -1';
                  lbl = '';
                  lbr = '';
              end
              extra.label{nextra} = ['Right Monoc' lbr];
              nextra = nextra+1;
              extra.label{nextra} = ['Left Monoc' lbl];
              if (isfield(Expt.Stimvals,'e3') & strcmp(Expt.Stimvals.e3,'ce')) | ~isempty(strfind(Expt.Header.Name,'ACNC')) | ~isempty(strfind(Expt.Header.Name,'DPxSF'))
                  splitextra(nextra-1:nextra) = 1;
              else
                  splitextra(nextra-1:nextra) = 0;
              end
          end
      end
  elseif strncmpi(str,'Uncorr',4)
      if isfield(Expt.Trials,'ce')
          nextra = nextra+1;
          extraexp{nextra} = 'ce] == 0';
          extra.label{nextra} = 'Uncorr';
          if (isfield(Expt.Stimvals,'e3') &  strcmp(Expt.Stimvals.e3,'ce')) | ~isempty(strfind(Expt.Header.Name,'ACNC')) | strcmp(Expt.Stimvals.e2,'c1')
              splitextra(nextra) = 1;
          else
              splitextra(nextra) = 0;
          end
      end
  elseif strmatch(str,'RandPhase')
    nextra = nextra+1;
    splitextra(nextra) = 0;
    extraexp{nextra} = 'rndphase] == 1';
    extra.label{nextra} = 'Random Phase';
  elseif strncmpi(str,'ShowN',5)
      showcounts = 1;
      psfargs = {psfargs{:}, 'showN'};
  elseif strncmpi(str,'xval',4)
      j = j+1;
      xvs = varargin{j};
      forcexvs = 1;
  elseif strncmpi(str,'yval',4)
      j = j+1;
      yval = varargin{j};
      ids = find([Expt.Trials.(btype)] == yval);
      if isempty(ids)
          fprintf('No Trials with %s =%.4f\n',btype,yval);
      else
          Expt.Trials = Expt.Trials(ids);
      end
  elseif strncmpi(str,'sval',4)
      j = j+1;
      ztype = varargin{j};
      j = j+1;
      zval = varargin{j};
      ids = find([Expt.Trials.(ztype)] == zval);
      if isempty(ids)
          fprintf('No Trials with %s =%.4f\n',btype,yval);
      else
          Expt.Trials = Expt.Trials(ids);
      end
  end
  j = j+1;
 end


if isfield(result,'ids')
    for j = 1:length(Expt.Trials)
        spktimes{j} = Expt.Trials(j).Spikes;
        counts(j) = sum(Expt.Trials(j).Spikes > 500 & Expt.Trials(j).Spikes < result.duration + 500);
    end
    if ~isfield(result.extras,'id')
        result.extras.id = [];
    end
    nc = 1;
    for k= 1:size(result.ids,2)
        for j= 1:size(result.ids,1)
            idx = result.ids{j,k};
            result.means(j,k) = mean(counts(idx));
            result.sd(j,k) = std(counts(idx));
            if gettimes
                result.spktimes{j,k} = cat(1, spktimes{idx});
            end
        end
        
        er = errorbar(result.x(:,k,nc),result.means(:,k,nc),result.sd(:,k,nc) ./ sqrt(max(result.n(:,k,nc),1)),'o-');
        set(er,'color',colors{k});
        hold on;
        plot(result.x(:,k),result.y(:,k),'-','color',colors{k},'marker','o','markerfacecolor',colors{k});
    end
    for k= 1:size(result.extras.id,2)
        idx = result.extras.id{k};
        result.extras.means(k) = mean(counts(idx));
        result.extras.sd(k) = std(counts(idx));
        plot(result.extras.x(k),result.extras.means(k),'o','color',colors{j});
        errorbar(result.extras.x(k),result.extras.means(k),result.extras.sd(k),'o','color',colors{j});
    end
end
if plotlfp
    fidx = find(ftfrq < lfprange(3));
    fidx = fidx(2:end);
    gidx = find(ftfrq < lfprange(2) & ftfrq > lfprange(1));
    agidx = find(ftfrq < lfprange(1) & ftfrq > lfprange(4));
    bgidx = find(ftfrq < lfprange(4) & ftfrq > lfprange(5));
    cgidx = find(ftfrq < lfprange(2));
end

Expt = FillExpt(Expt, 'st');
Expt = FillExpt(Expt, type); %make sure stimval is set in case n == 1
stimtype = GetEval(Expt,'st','mode');
[eye, eyes] = GetEval(Expt,'me','mode');
if isempty(duration)
    duration = prctile([Expt.Trials.End] - [Expt.Trials.Start],50);
else
    id = find([Expt.Trials.End]- [Expt.Trials.Start] > duration * 0.98);
    Expt.Trials = Expt.Trials(id);
end

if isfield(Expt.Trials,'RespDir') & mean(abs([Expt.Trials.RespDir])) > 0.5 & includebad ~= 1
    id = find([Expt.Trials.RespDir] ~= 0);
    Expt.Trials = Expt.Trials(id);
end

freq = GetEval(Expt,'tf')/10000;
if(freq > 0)
    period = 1/freq;
else
    period = duration;
end

if mksdf
    if isempty(times)
        if mksdf == 2
            times = -1000:20:4000;
        else
            times = -preperiod:20:duration+postperiod;
        end
    end
    result.sdftimes = times;
elseif getmod
    times = 500:1:duration+500;
end

result.type{1} = type;
if exist('btype','var')
    result.type{2} = btype;
else
    result.type{2} = 'e0';
end
if exist('ctype','var')
result.type{3} = ctype;
end

mucounts = [];
for j = 1:length(Expt.Trials);
    counts(j) = length(find([Expt.Trials(j).Spikes] > latency & ...
        [Expt.Trials(j).Spikes] < duration+latency));
    if gettimes
        spktimes{j} = Expt.Trials(j).Spikes;
    end
    
    if isfield(Expt.Trials,'OSpikes')
        mucounts(j) = length(find([Expt.Trials(j).OSpikes] > latency & ...
            [Expt.Trials(j).OSpikes] < duration+latency));
    end
    Expt.Trials(j).count = counts(j);
end

if showmu == 0
    bcounts = mucounts;
else
    bcounts = [];
end

if acov
    Expt = FillTrials(Expt,'tf');
    [ac, acoveach] = autocorrelate(Expt.Trials,latency, duration);
else
    acoveach = [];
end
if dosqrt
    counts = counts .^0.5;
end
if showrates
    counts = counts * 10000/duration;
end

bid = find([Expt.Trials.st] == 0);
if length(bid) && length(bid) < length(Expt.Trials)/3  %blank interleaved
    [Expt.Trials(bid).(type)] = deal(-1009);
end
if isfield(Expt.Trials,'ce')
    bid = find([Expt.Trials.ce] == 0);
    if length(bid) && length(bid) < length(Expt.Trials)/3  %uncorr interleaved
        [Expt.Trials(bid).(type)] = deal(-1005);
    end
end
if ~forcexvs 
    if isempty(xvs)
        xvs = sort(unique([Expt.Trials.(type)]));
        xvs = xvs(find(~isnan(xvs)));
        if strcmp(type,'dx') && max(xvs) > 100  %%mistaken highx interleave
            id = find([Expt.Trials.(type)] >= 100);
            [Expt.Trials(id).(type)] = deal(NaN);
            xvs = xvs(find(xvs < 100));
        end
    end
    if logx
        mindiff = 0;
    else
        [xvs, mindiff] = RemoveDuplicates(xvs,Expt);
    end
else
    mindiff = 0;
end
nn = 0;
nsum = 0;
%First cacluate all the interleaved extras. Keep a tally of all
%these trials so that they can be excluded from the main expt. 
%start off with a list of trials to be excluded altogether
extraidx = [];

extraidx = find([Expt.Trials.Trial] < 0);
result.excluded = extraidx;
for ie = 1:length(extraexp)
    extra.n(ie) = 0;
    if ~isempty(extraexp{ie})
%first identify trials that meet condition 
        idx = eval(['find([Expt.Trials.' extraexp{ie} ');']);
%then exclude any that have already been included in earlier condition
%(makes sure a blank stim is not also counted as uncorr, for example)
        idx = setdiff(idx,extraidx);
% if most of all trials are and "extra", it must be a real stim condition
% (e.g. no stim for long spontaneous measures)
        if length(idx) >= nmin && length(idx) < length(Expt.Trials) * 0.9
            extra.x(ie) = mean(xvs(find(xvs > -1000)));
            extra.means(ie) = mean(counts(idx));
            extra.n(ie) = length(counts(idx));
            extra.sd(ie) = std(counts(idx));
            extra.id{ie} = idx;
            extra.ids{ie} = [Expt.Trials(idx).id];
            extra.counts{ie} = counts(idx);
            extraidx = [extraidx idx];
            if gettimes
                extra.spktimes{ie} = cat(1,spktimes{idx});
            end
            if ~isempty(acoveach)
                extra.acov{ie} = mean(acoveach(idx,:));
            end
            if ~isempty(idx)
                nsum = nsum + length(idx);
                nn = nn + 1;
            end
        end
    end
    if extra.n(ie) == 0 
        extra.means(ie) = NaN;
        extra.x(ie) = mean(xvs);
    end
end


if plotlfp  && isfield(Expt.Trials,'FTlfp')
    if extra.n(1) > 0
        idx = extra.id{1};
        fts = abs([Expt.Trials(idx).FTlfp]);
        ft = mean(fts,2);
        aft = mean(fts,1);
        if length(ft) >= max(agidx);
        blanklfpwr = mean(ft(gidx));
        ablanklfpwr = mean(ft(agidx));
        bblanklfpwr = mean(ft(bgidx));
        blanklfpse = std(mean(fts(gidx,:),1))/sqrt(length(idx));
        ablanklfpse = std(mean(fts(agidx,:),1))/sqrt(length(idx));
        bblanklfpse = std(mean(fts(bgidx,:),1))/sqrt(length(idx));
%Calculate mean resp of all non-blank stimuli.
        aidx = setdiff([1:length(Expt.Trials)],extra.id{1});
        oft = mean(abs([Expt.Trials(aidx).FTlfp]),2);
        respf = smooth(oft(fidx)-ft(fidx),3);
        result.lfpower = smooth(oft,5);
        [a, peakf] = max(respf);
        j = peakf;
        while(respf(j) > a/5 & j < length(ftfrq))
            j = j+1;
        end
        maxf = j-1;
        j = peakf;
        while(respf(j) > a/5 & j > 1)
            j = j-1;
        end
        minf = j+1;
        cgidx = minf:maxf;
        cblanklfpwr = mean(ft(cgidx));
        lfpautof = [ftfrq(minf) ftfrq(maxf)];
        end
    else
        blanklfpwr = 0;
        ablanklfpwr = 0;
        bblanklfpwr = 0;
        cblanklfpwr = 0;
        lfpautof = [];
        oft = mean(abs([Expt.Trials.FTlfp]),2);
        result.lfpower = smooth(oft,5);
    end
end


idx = setdiff(1:length(Expt.Trials),extraidx);

if ~forcexvs && ~logx
    xvs = sort(unique([Expt.Trials(idx).(type)]));
    xvs = xvs(find(~isnan(xvs)));
    [xvs, mindiff] = RemoveDuplicates(xvs, Expt);
end
if exist('btype','var')
  bvals = eval(['sort(unique([Expt.Trials(idx).' btype ']));']);
else
  bvals = NaN;
end
if GetEval(Expt,'st') == 10 & isfield(Expt.Trials,'st')
    id = find([Expt.Trials.st] == 3)
    sfs = unique([Expt.Trials(id).sf]);
    
    ida = find([Expt.Trials(id).sf] == sfs(1));
    idb = find([Expt.Trials(id).sf] == sfs(2));
    bvals = [bvals -1001];
    [Expt.Trials(id(ida)).(btype)] = deal(-1001);
    bvals = [bvals -1002];
    [Expt.Trials(id(idb)).(btype)] = deal(-1002);
end

if exist('ctype','var')
  cvals = eval(['sort(unique([Expt.Trials(idx).' ctype ']));']);
  if strmatch(ctype,'ce')
      cvals = fliplr(cvals);
  end
else
  cvals = NaN;
end

result.nlines = length(bvals);
result.duration = duration;
h = [];

if autonmin
    for nc = 1:length(cvals)
        if ~isnan(cvals(nc))
            cidx  = eval(['find([Expt.Trials.' ctype '] == cvals(nc));']);
        else
            cidx = 1:length(Expt.Trials);
        end
    for ie = 1:length(bvals)
        if ~isnan(bvals(ie))
            bidx  = eval(['find([Expt.Trials.' btype '] == bvals(ie));']);
        else
            bidx = 1:length(Expt.Trials);
        end
        for x = xvs;
            idx = find(abs([Expt.Trials.(type)]- x) < mindiff & [Expt.Trials.st] > 0);
            idx = intersect(idx,bidx);
            idx = intersect(idx,cidx);
            idx = setdiff(idx,extraidx);
            if ~isempty(idx)
                nsum = nsum + length(idx);
                nn = nn + 1;
            end
        end
    end
end
    autonmin = nsum/nn;
    nmin = ceil(autonmin/2);
end

tic;
step = 0.1; %%default offset
if logx
    stepratio = (max(xvs)./min(xvs))^0.025;
end

result.allhandles = [];

for nc = 1:length(cvals)
    if ~isnan(cvals(nc))
        cidx  = eval(['find([Expt.Trials.' ctype '] == cvals(nc));']);
    else
        cidx = 1:length(Expt.Trials);
    end
    if exist('btype','var')
    for ie = 1:length(bvals)
        if ~isnan(bvals(ie))
            bidx  = eval(['find([Expt.Trials.' btype '] == bvals(ie));']);
            nbval(ie) = length(bidx);
        else
            nbval(ie) = sum(isnan([Expt.trials.(btype)]));
        end
    end
    bid = find(nbval < mean(nbval/5)); %%do't allow conditions with many fewer repeats
    bid = find(nbval > mean(nbval/5)); %%do't allow conditions with many fewer repeats
    bvals = bvals(bid);
    end
    for ie = 1:length(bvals);
        ix = 1;
        if ~isnan(bvals(ie))
            bidx  = eval(['find([Expt.Trials.' btype '] == bvals(ie));']);
        else
            bidx = 1:length(Expt.Trials);
        end
        result.linevals(ie,nc) = SetLineval(bvals,bvals(ie));
        result.colors{ie,nc} = colors{ie+addn};
        step = (max(xvs)-min(xvs))/50;

        % combined files might have .st in one and not antoher
        if isfield(Expt.Trials,'st') & length([Expt.Trials.(type)]) ~= length([Expt.Trials.st])
            for j = 1:length(Expt.Trials)
                if isempty(Expt.Trials(j).st)
                    Expt.Trials(j).st = Expt.Stimvals.st;
                end
            end
        end

        newcode = 1;
        ccidx = setdiff(intersect(bidx,cidx),extraidx);
        for xi = 1:length(xvs)
            x = xvs(xi);
            if newcode
                idx = find(abs([Expt.Trials(ccidx).(type)]- x) <= mindiff);
                idx = ccidx(idx);
            else
                idx = find(abs([Expt.Trials.(type)]- x) <= mindiff & [Expt.Trials.st] > 0);
                if Expt.Stimvals.st ~= 0 %should be covered by extraidx
                    idx = setdiff(idx,[Expt.Trials.st] == 0);
                end
                idx = find([Expt.Trials.(type)] == x & [Expt.Trials.st] > 0);
                idx = cidx(idx);
                idx = intersect(idx,bidx);
                idx = intersect(idx,cidx);
                idx = setdiff(idx,extraidx);
            end
            if (length(idx) >= nmin | fillall) & (x > 0 | ~logx)
                if x < -999
                    x = SetLineval(xvs,x);
                end
                allcounts{ix,ie,nc} = counts(idx);
                result.x(ix,ie,nc) = x;
                result.y(ix,ie,nc) = result.linevals(ie);
                result.z(ix,ie,nc) = cvals(nc);
                result.n(ix,ie,nc) = length(counts(idx));
                if isfield(Expt.Trials,'thinor')
                    result.thinor(ix,ie,nc) = mean([Expt.Trials(idx).thinor]);
                end
                if isempty(idx)
                    result.means(ix,ie,nc) = NaN;
                    result.sd(ix,ie,nc) = NaN;
                else
                    if ~isempty(bcounts)
                        result.bmeans(ix,ie,nc) = mean(bcounts(idx));
                        result.bsd(ix,ie,nc) = std(bcounts(idx));
                    end
                    result.means(ix,ie,nc) = mean(counts(idx));
                    result.sd(ix,ie,nc) = std(counts(idx));
                end
                result.flags(ix,ie,nc) = 0;
                result.counts{ix,ie,nc} = counts(idx);
                result.ids{ix,ie,nc} = [Expt.Trials(idx).id];
                result.tidx{ix,ie,nc} = idx;
                if gettimes
                    result.spktimes{ix,ie,nc} = cat(1, spktimes{idx});
                end
    if isfield(Expt.Trials,'RespDir')
                    result.RespDir{ix,ie,nc} = [Expt.Trials(idx).RespDir];
                end
                if ~isempty(acoveach)
                    ac = mean(acoveach(idx,:));
                    result.acov(ix,ie,:) = ac./ac(1);
                    result.acov(ix,ie,1) = 0;
                    result.period(ix,ie) = 10000./mean([Expt.Trials(idx).tf]);
                end
                if isfield(Expt.Trials,'mf')
                    freq = mean([Expt.Trials(idx).mf])/10000;
                    period = 1/freq;
                elseif isfield(Expt.Trials,'tf') & length(idx) > 0
                    freq = mean([Expt.Trials(idx).tf])/10000;
                    if freq > 0
                        period = 1/freq;
                    else
                        period = duration;
                    end
                end
                if periodic
                    sdfargs{1} = 'period';
                    sdfargs{2} = period;
                    sdfw = period/periodsdf;
                end
                if mksdf == 1
                    result.sdfs{ix,ie,nc} = trigsdfa(Expt.Trials(idx),sdfw,times,sdfargs{:});
                    if exist('btype','var')
                        result.labels{ix,ie,nc} = sprintf('%s %f %s %f',btype,bvals(ie),type,x);
                    else
                        result.labels{ix,ie,nc} = sprintf('%s %.3f(%d)',type,x,result.n(ix,ie,nc));
                    end
                elseif mksdf == 2
                    [result.sdfs{ix,ie,nc} result.nsacs(ix,ie,nc)] = mksacsdf(Expt.Trials(idx),sdfw,times,'halfgauss',sdfargs{:});
                    if exist('btype','var')
                        result.labels{ix,ie,nc} = sprintf('%s %f %s %f',btype,bvals(ie),type,x);
                    else
                        result.labels{ix,ie,nc} = sprintf('%s %.3f(%d)',type,x,result.nsacs(ix,ie,nc));
                    end
                elseif getmod & ~isempty(idx)
                    sdf = trigsdfa(Expt.Trials(idx),sdfw, ...
                        times,'raw',sdfargs{:});
                    result.sdfs{ix,ie,nc} = sdf;
                    % * 1000 to convert spikes/0.1ms to spikes/sec;
                    result.f1(ix,ie,nc) = (10000/length(idx)) * famp(times(1:length(sdf)),sdf,freq);
                end

                lfpch = 1;
                if plotlfp
                    if length(ft) == 0
                        plotlfp = 2;
                    end
                    subplot(2,1,2);
                    fts = abs([Expt.Trials(idx).FTlfp]);
                    ft = mean(abs([Expt.Trials(idx).FTlfp]),2);
                    if exist('btype','var')
                        stimlab = [val2str(result.x(ix,ie,nc),type,stimtype,Expt,[]) val2str(result.y(ix,ie,nc),btype,stimtype,Expt,[])];
                    else
                        stimlab = [val2str(result.x(ix,ie,nc),type,stimtype,Expt,[])];
                    end
                    if isempty(idx)
                        fprintf('No Data for %s\n',stimlab);
                    else
                        if plotlfp == 2 %time domain
                            addl = mod(addl,100);
                            ts = [1:size([Expt.Trials.LFP],1)] .* Expt.Header.LFPsamplerate;
                            lfph(nlfp) = plot(ts,mean([Expt.Trials(idx).LFP],2),'color',colors{ix+addn+addl});
                            result.lfpt(ix,ie,:) = mean([Expt.Trials(idx).LFP],2);
                        else
                            lfph(nlfp) = plot(ftfrq(fidx),ft(fidx),'color',colors{ix+addn+addl});
                        end
                        lfplabels{nlfp} = stimlab;
                        nlfp = nlfp + 1;
                        if length(ft) > 1
                        result.lfpwr(ix,ie,nc) = mean(ft(gidx)) - blanklfpwr;
                        result.alfpwr(ix,ie,nc) = mean(ft(agidx)) - ablanklfpwr;
                        result.blfpwr(ix,ie,nc) = mean(ft(bgidx)) - bblanklfpwr;
                        result.clfpwr(ix,ie,nc) = mean(ft(cgidx)) - cblanklfpwr;
                        result.lfpwrs{ix,ie,nc} = mean(fts(gidx,:),1) - blanklfpwr;
                        result.alfpwrs{ix,ie,nc} = mean(fts(agidx,:),1) - ablanklfpwr;
                        result.blfpwrs{ix,ie,nc} = mean(fts(bgidx,:),1) - bblanklfpwr;
                        result.clfpwrs{ix,ie,nc} = mean(fts(cgidx,:),1) - cblanklfpwr;
                        result.lpfsnr = Expt.Header.LFPsnr;
                        end
                        hold on;
                    end
                    if(laxis ~= 0)
                        axes(laxis)
                    else
                        subplot(2,1,1);
                        hold off;
                    end
                end
                if triglfp
                    lw = 200;
                    if Expt.Header.lfplen < 401
                        lw = (length(Expt.Trials(1).LFP)/2)-1;
                    else
                        lw = 200;
                    end
                    [result.triglfp(ix,:), nspk] = SpTrigLFP(Expt.Trials(idx),duration, Expt.Header.CRsamplerate, lw);
                    result.trignspk(ix) = nspk;
                end
                ix = ix+1;
            end %% end of if isempty(idx)
            addl = addl+length(bvals);
        end




        if(laxis ~= 0)
            axes(laxis)
        end

% This test needs && so that no field x stops testing
        if showplot && ~mksdf && isfield(result,'x') && ~plotsimple && ie <= size(result.x,2) && nc <= size(result.x,3)
            result.means(find(result.n == 0)) = NaN;
            is = mod(ie-1,length(symbols))+1;
            if isempty(colorids)
                c = colors{ie+addn};
            else
                c = colors{colorids(ie,nc)};
            end
            if showplot == 2
                allh = [];
                for k = 1:length(result.x(:,ie,nc))
                    ph = plot(result.x(k,ie,nc),result.counts{k,ie,nc},symbols(is), ...
                        'color',colors{ie+addn});
                    allh = [allh; ph];
                    hold on;
                    h(ie) = allh(1);
                end
            else
                h(ie) = plot(result.x(:,ie,nc),result.means(:,ie,nc),symbols(is), ...
                    'color',c);
                result.colors{ie,nc} = c;
                allh = h(ie);
                result.handles(ie) = h(ie);
                result.allhandles = [result.allhandles h(ie)];
                if isfield(result,'bmeans')
                    hold on;
                    bh(ie) = plot(result.x(:,ie,nc),result.bmeans(:,ie,nc),symbols(is), ...
                        'color',c);
                end
            end

            if fillsymbols & nc == 1
                set(allh,'MarkerFaceColor',c);
            end
            hold on;
            if getmod
                plot(result.x(:,ie,nc),result.f1(:,ie,nc),'o-', ...
                    'color',colors{ie+addn});
            end
            er = errorbar(result.x(:,ie,nc),result.means(:,ie,nc),result.sd(:,ie,nc) ./ sqrt(max(result.n(:,ie,nc),1)),'o-');
            set(er,'color',c,'linestyle',linestyle{nc});
            result.allhandles = [result.allhandles er];
            if nolines
                set(er,'linestyle','none');
            end
            if isfield(result,'bmeans')
                ber = errorbar(result.x(:,ie,nc),result.bmeans(:,ie,nc),result.bsd(:,ie,nc) ./ sqrt(max(result.n(:,ie,nc),1)),'o-');
                set(ber,'color',colors{ie+addn},'linestyle',linestyle{nc});
            end

            if(plotlfp ==1)
                for li = 1:size(result.lfpwrs,1)
                    lstd(li) = std(result.lfpwrs{li,ie,nc}) ./ sqrt(length(result.lfpwrs{li,ie,nc}));
                    astd(li) = std(result.alfpwrs{li,ie,nc}) ./ sqrt(length(result.alfpwrs{li,ie,nc}));
                    bstd(li) = std(result.blfpwrs{li,ie,nc}) ./ sqrt(length(result.blfpwrs{li,ie,nc}));
                    cstd(li) = std(result.blfpwrs{li,ie,nc}) ./ sqrt(length(result.blfpwrs{li,ie,nc}));
                end
                result.lfpautorange = lfpautof;
                laxis = gca;
                if(raxis == 0)
                    raxis = axes('Position',get(gca,'Position'),'color','none');
                    hold off;
                else
                    axes(raxis);
                end
                errorbar(result.x(:,ie,nc),result.lfpwr(:,ie,nc),lstd);
                hold on;
                plot(result.x(:,ie,nc),result.lfpwr(:,ie,nc), 'o-','color',colors{ie+addn});
                errorbar(result.x(:,ie,nc),result.alfpwr(:,ie,nc),astd);
                plot(result.x(:,ie,nc),result.alfpwr(:,ie,nc),'^-','color',colors{ie+addn});
                errorbar(result.x(:,ie,nc),result.blfpwr(:,ie,nc),bstd);
                plot(result.x(:,ie,nc),result.blfpwr(:,ie,nc), 's-','color',colors{ie+addn});
                errorbar(result.x(:,ie,nc),result.clfpwr(:,ie,nc),cstd);
                plot(result.x(:,ie,nc),result.clfpwr(:,ie,nc), '*-','color',colors{ie+addn});
                set(raxis,'Color','None','Xlim',get(gca,'Xlim'),'YAxisLocation','Right','Xscale',get(laxis,'Xscale'));
                axes(laxis);
            end
        else
            result.handles(1) = 0;
        end



        if isnan(bvals(ie))
            if strmatch(type,'cvsign')
                labels{ie} = sprintf('Or%s',sprintf(' %.0f',unique([Expt.Trials.or]))); 
            else
                labels{ie} = val2str(eye,'me',stimtype,Expt,0);
            end
        else
            if range(bvals) < 0.1
                nfp = 3;
            else
                nfp = 2;
            end
            labels{ie} = val2str(bvals(ie),btype,stimtype,Expt,nfp);
        end
        if showplot & result.handles(1)
            hold on;
        end

    end
    if exist('ctype','var')
        result.(ctype)(nc) = cvals(nc);
    end
end %end of nc in cvals loop

result.x2 = 1;
if isfield(Expt.Stimvals,'x2') && Expt.Stimvals.x2 == 0 && Expt.Stimvals.n2 > 0
    if sethold == 0
        hold off;
    end
    result.x2 = 0;
   for nc = 1:size(result.x,3)
       nx = sum(result.n(:,:,nc),1);
       ny = sum(result.n(:,:,nc),2);
       [a,xi] = max(nx);
       [a,yi] = max(ny);
       er = errorbar(result.x(:,xi,nc),result.means(:,xi,nc),result.sd(:,xi,nc) ./ sqrt(max(result.n(:,xi,nc),1)),'o-');
       hold on;
       set(er,'color',colors{1},'linestyle',linestyle{nc});
       allh = plot(result.x(:,xi,nc),result.means(:,xi,nc),symbols(nc));
       result.allhandles = [result.allhandles er];
       if fillsymbols & nc == 1
           set(allh,'MarkerFaceColor',colors{1});
       end
       
       
       er = errorbar(result.y(yi,:,nc),result.means(yi,:,nc),result.sd(yi,:,nc) ./ sqrt(max(result.n(yi,:,nc),1)),'o-');
       set(er,'color',colors{2},'linestyle',linestyle{nc});
       result.allhandles = [result.allhandles er];
       allh = plot(result.y(yi,:,nc),result.means(yi,:,nc),symbols(nc));
       if fillsymbols & nc == 1
           set(allh,'MarkerFaceColor',colors{2});
       end
       
   end
end

if isfield(result, 'n') %% empty if all stims are "extra"
    id = find(result.n == 0);
    result.means(id) = NaN;
end


%Extra interleaved stimul. extraexp defines the 
if isempty(ie)
    ie = 0;
end
nlin = ie;
nplots = ie+1;
nsplit = 1;
extraoffset = [1 -1 2 -2];
for ie = 1:length(extraexp)
if extra.n(ie) > 0 & ~plotsimple
    if splitextra(ie) & showplot
        delta = extraoffset(nsplit) * mean(diff(xvs))/5;
        for k = 1:length(bvals)
            idx = find([Expt.Trials.(btype)] == bvals(k));
            idx = intersect(idx,extra.id{ie});
            resp = mean(counts(idx));
            extra.splitmeans(ie,k) = resp;
            allh(k) = plot(extra.x(ie)+delta,resp,'o', ...
                'color',colors{k+addn});
            er = errorbar(extra.x(ie)+delta,resp,std(counts(idx)) ./ sqrt(length(idx)),'o-');
            set(er,'color',colors{k+addn});
            if(showcounts & showplot)
                if logx
                text(extra.x(ie)+step+delta,resp,sprintf('%d',length(idx)),'color',colors{k+addn});
                else
                text(extra.x(ie)+step+delta,resp,sprintf('%d',length(idx)),'color',colors{k+addn});
                end
            end
        end
        nsplit = nsplit+1;
        h(nplots) = allh(1);
    elseif showplot == 2
        allh = plot(extra.x(ie),extra.counts{ie},'o', ...
            'color',colors{ie+nlin+addn});
        h(nplots) = allh(1);
        er = errorbar(extra.x(ie),extra.means(ie),extra.sd(ie) ./ sqrt(extra.n(ie)),'o-');
        set(er,'color',colors{ie+nlin+addn});
        if(showcounts & showplot)
            text(extra.x(ie)+step,extra.means(ie),sprintf('%d',extra.n(ie)),'color',colors{ie+nlin});
        end
    elseif showplot
        cid = 1+mod(ie+nlin+addn-1,length(colors));
        h(nplots) = plot(extra.x(ie),extra.means(ie),'o', ...
            'color',colors{cid});
        er = errorbar(extra.x(ie),extra.means(ie),extra.sd(ie) ./ sqrt(extra.n(ie)),'o-');
        set(er,'color',colors{cid});
        if(showcounts & showplot)
            text(extra.x(ie)+step,extra.means(ie),sprintf('%d',extra.n(ie)),'color',colors{cid});
        end
    end

  labels{nplots} = extra.label{ie};
  
  if plotlfp          
      idx = extra.id{ie};
      if showplot
          subplot(2,1,2);
      end
      fts = abs([Expt.Trials(idx).FTlfp]);
      ft = mean(abs([Expt.Trials(idx).FTlfp]),2);
      if length(ft) > 1
      extra.lfpwr(ie) = mean(ft(gidx)) - blanklfpwr;
      extra.lfpwrs{ie} = mean(fts(gidx,:),1) - blanklfpwr;
      extra.alfpwrs{ie} = mean(fts(agidx,:),1) - ablanklfpwr;
      extra.blfpwrs{ie} = mean(fts(bgidx,:),1) - bblanklfpwr;
      extra.clfpwrs{ie} = mean(fts(cgidx,:),1) - cblanklfpwr;
      extra.alfpwr(ie) = mean(ft(agidx)) - ablanklfpwr;
      extra.blfpwr(ie) = mean(ft(bgidx)) - bblanklfpwr;
      extra.clfpwr(ie) = mean(ft(cgidx)) - cblanklfpwr;
      extra.lfpse(ie) = std(mean(fts(gidx,:),1))/sqrt(length(idx));
      extra.alfpse(ie) = std(mean(fts(agidx,:),1))/sqrt(length(idx));
      extra.blfpse(ie) = std(mean(fts(bgidx,:),1))/sqrt(length(idx));
      extra.clfpse(ie) = std(mean(fts(cgidx,:),1))/sqrt(length(idx));
      end
      if showplot & plotlfp == 1
          if ie ==1  %%Blank stimulus
              plot(ftfrq(fidx),ft(fidx),'color','k','linewidth',2);
              idx = setdiff([1:length(Expt.Trials)],extra.id{ie});
              oft = mean(abs([Expt.Trials(idx).FTlfp]),2);
              bar(ftfrq(fidx),smooth(oft(fidx)-ft(fidx),3));
              axes(raxis);
%Draw a dashed horizontal line for the LFP response to Blank.  Add a little
%to Xlim to avoid plotting at x = 0 in case axes turn logrithmic later.
          
              plot(get(gca,'Xlim')+0.001,[extra.lfpwr(ie) extra.lfpwr(ie)],':','color','k');
              subplot(2,1,2);
          else
              plot(ftfrq(fidx),ft(fidx),'color',colors{ix+addn});
          end
          hold on;
          if(laxis ~= 0)
              axes(laxis)
          else
              subplot(2,1,1);
          end
      elseif plotlfp == 2
          ts = [1:size([Expt.Trials.LFP],1)] .* Expt.Header.LFPsamplerate;
          plot(ts,mean([Expt.Trials(idx).LFP],2),'color','k');
      end
  end
  nplots = nplots+1;
end
end

if isfield(result,'n') && result.x2 > 0%% empty if no data
id = find(sum(result.n) > 0);
if length(id) < size(result.x,2) & ~fillall
result.x = result.x(:,id);
result.y = result.y(:,id);
result.n = result.n(:,id);
result.means = result.means(:,id);
result.ids = result.ids(:,id);
result.counts = result.counts(:,id);
result.sd = result.sd(:,id);
end
end


if isfield(result,'x') && result.x2 > 0
%check rows have consistent X vals. Can have 0s if missing data....
for j = 1:size(result.x,1)
    if length(unique(result.x(j,:))) > 1
        id = find(result.n(j,:) > 0);
        xs = unique(result.x(j,id));
        if length(xs) == 1
            result.x(j,:) = x;
        end
    end
end
end


if plotlfp == 1 & showplot
    subplot(2,1,2);
    if isempty(result.lfpautorange)
        title(sprintf('No Blank'));
    else
        title(sprintf('%.0f - %.0fHz (SNR %.1f)',result.lfpautorange(1),result.lfpautorange(2),Expt.Header.LFPsnr));
    end        
    if nlfp > 1 & nlfp < 10        
        legend(lfph,lfplabels);
    end
end

if mksdf
    if sethold == 0
    hold off;
    end
    n = 1;
    cn = 0; %to allow cycling through colors in case too many lines
    if isempty(subsdf)
    showsdf = ones(size(result.sdfs));
    else
    showsdf = zeros(size(result.sdfs));
    id = eval(['find(' subsdf ')']);
    showsdf(id) = 1;
    end
    for j = 1:size(result.sdfs,1)
    for k = 1:size(result.sdfs,2)
        if ~isempty(result.sdfs{j,k}) & showsdf(j,k)
            if periodic
                h(n) = plot([1:length(result.sdfs{j,k})]./length(result.sdfs{j,k}),result.sdfs{j,k},'color',colors{n});
            elseif(length(times) == length(result.sdfs{j,k}))
                h(n) = plot(times,result.sdfs{j,k},'color',colors{n-cn},'linestyle',linestyle{nline});
            else
                h(n) = plot(result.sdfs{j,k},'color',colors{n-cn});
            end
            if isfield(Expt.Trials,'RespDir')
                labels{n} = sprintf('%s RespDir %.1f(%.2f)',result.labels{j,k},mean([Expt.Trials.RespDir]),mean([Expt.Trials.(type)]));
            else
            labels{n} = result.labels{j,k};
            end
          n = n + 1;
          if n > length(colors)
              cn = n-1;
          end
          hold on;
        end
      end
  end
  legend
  if showsum
      result.sdfall = trigsdfa(Expt.Trials,sdfw,times,sdfargs{:});
      if length(times) == length(result.sdfall)
          h(n) = plot(times,result.sdfall,'k','linewidth',2);
      else
          h(n) = plot(result.sdfall,'k','linewidth',1);
      end
      labels{n} = 'All';
  end
  if isfield(Expt.Trials,'dur')
  endtime = median([Expt.Trials.dur]);
  plot([endtime endtime],get(gca,'ylim'),'k:');
  end
end

if showplot
if isfield(Expt,'Names') & isfield(Expt.Names,'Code')
    id = strmatch(type,{Expt.Names.Code});
    if ~isempty(id)
        xlabel(Expt.Names(id).Label);
    else
        xlabel(type);
    end
else
    xlabel(type);
end

%Make Sure axis label is on correct one if plotting LFP/
if laxis
    axes(laxis);
end
if showrates
    ylabel(sprintf('Spike Rate'));
else
    ylabel(sprintf('Spike Count (%.2f sec)',duration/10000));
end


if isfield(result,'acov')
    if plotsimple
        ylabel('P(spike)');
        xlabel('t (cycles)');
        nplots = 1;
    end
    off = round(result.period./4);
    for j = 1:size(result.period,1)
        for k = 1:size(result.period,2)
            ts = [1:result.period(j,k)]+off(j,k);
            [a,b] = famp(ts,squeeze(result.acov(j,k,ts)),1/result.period(j,k));
            Fn(j,k,1) = real(b); 
            [a,b] = famp(ts,squeeze(result.acov(j,k,ts)),2/result.period(j,k));
            Fn(j,k,2) = real(b); 
            Fn(j,k,3) = mean(result.acov(j,k,:));
            
            if plotsimple & showplot
                sm = smooth(squeeze(result.acov(j,k,1:max(ts)*2)),off(j,k)/4,'gauss');
                h(nplots) = plot([1:max(ts)*2]./result.period(j,k),sm,'color',colors{j*k});
                labels{nplots} = sprintf('TF %.2f',10000./result.period(j,k));
                nplots = nplots+1;
                hold on;
            end
        end
    end
% put F2 on the demoninatior - this is usually nonzero even for simple
% cells. Ratio > 1 = simple cell
    result.f1f2 = Fn(:,:,1)./Fn(:,:,2);
    result.f1f2sum = (Fn(:,:,1)+Fn(:,:,2))./Fn(:,:,3);
    [a,b] = max(result.f1f2sum(:)); % strongest modulation
    result.simple = [result.f1f2(b) result.f1f2sum(b)];
    if plotsimple
        if result.simple(1) > 1
            ctype = 'Simple';
        else
            ctype = 'Complex';
        end
        result.title = title(sprintf('%s %s (F1/F2 = %.2f, Modulation %.2f)',...
            splitpath(Expt.Header.Name),ctype,result.simple(1),result.simple(2)));
        set(gca,'xscale','lin','yscale','lin');
        logx = 0;
        result.handles = h;
    end 
%    plot(Fn(:,:,2)+Fn(:,:,1),Fn(:,:,2)./Fn(:,:,1),'o');
end

end
result.toc = toc;
    
result.legend.h = h;
result.legend.labels = labels;
if((isstr(legendpos) | legendpos < 6) & showplot)
    hid = find(ishandle(h) & h > 0);
    if isstr(legendpos)
        legend(h(hid),labels{hid},'Location',legendpos);
    else
        legend(h(hid),labels{hid},legendpos);
    end
end



if(logx & ~mksdf)
  set(gca,'Xscale','log');
end


%Need to do this in case hold is on and Ymax needs to change. 
if showplot
    set(gca,'Ylimmode','auto');
    Yl = get(gca,'Ylim');
    Yl(1) = 0;
    set(gca,'Ylim',Yl);
    Expt.Plot.Ymax = Yl(2);
    if(raxis ~= 0)
        axes(raxis);
        set(raxis,'Color','None','Xlim',get(laxis,'Xlim'),'YAxisLocation','Right','Xscale',get(laxis,'Xscale'));
        ylabel(sprintf('LFP %d -%dHz',lfprange(1),lfprange(2)));
    end

    if exist('btype','var') &&  strcmp(btype,'xydir')
        yl = get(gca,'Ylim');
        xtick = get(gca,'xtick');
        yoff = Expt.Stimvals.rf(2)-Expt.Stimvals.rf(1);
        ytick = xtick + yoff;
        ytick = ceil(ytick .*20)./20;
        for j = 1:length(ytick)
            text(ytick(j)-yoff,yl(1)+diff(yl)/30,sprintf('%.2f',ytick(j)),'color','b','HorizontalAlignment','center');
        end
    end
end
[result.name, result.expname] = GetEval(Expt, 'name');
result.extras = extra;

if triglfp
    subplot(2,1,1);
    hold off;
    t = ([1:size(result.triglfp,2)]-size(result.triglfp,2)/2)./(10 * Expt.Header.CRsamplerate);
    for j = 1:size(result.triglfp,1)
        plot(t,result.triglfp(j,:),'color',colors{j});
        labels{j} = sprintf('%d',result.trignspk(j));
        hold on;
    end
    if isfield(Expt.Trials,'RespDir')
        ida = find([Expt.Trials.RespDir] == 1);
        [result.lfpchoice{1}, nspk] = SpTrigLFP(Expt.Trials(ida),duration, Expt.Header.CRsamplerate, lw);
        labels{j+1} = sprintf('%d',nspk);
        idb = find([Expt.Trials.RespDir] == -1);
        [result.lfpchoice{2}, nspk] = SpTrigLFP(Expt.Trials(idb),duration, Expt.Header.CRsamplerate, lw);
        labels{j+2} = sprintf('%d',nspk);
        plot(t,result.lfpchoice{1},'color',colors{j+1});
        plot(t,result.lfpchoice{2},'color',colors{j+2});
        cft = abs(fft(result.lfpchoice{1})) + abs(fft(result.lfpchoice{2}));
    else
        [result.lfpchoice{1}, nspk] = SpTrigLFP(Expt.Trials,duration, Expt.Header.CRsamplerate, lw);
        plot(t,result.lfpchoice{1},'color','k','linewidth',2);
        cft = abs(fft(result.lfpchoice{1}));
    end
    legend(labels);
    raxis = 0;
    subplot(2,1,2);
    plot(ftfrq(fidx),result.lfpower(fidx),'k','linewidth',2);
    cffrq = (0:length(cft)-1)/(length(cft)/(10000 * Expt.Header.CRsamplerate));
    cft = -cft .* max(result.lfpower(fidx))/max(cft);
    cfid = find(cffrq < ftfrq(fidx(end)));
    
    plot(cffrq(cfid),cft(cfid),'r','linewidth',2);
end
if showpcolor & size(result.x,2) > 1
    [X,Y,Z] = fillpmesh(result.x, result.y, result.means);
    step(1) = 0.5;
    for ie = 1:length(result.extras.means)
        if splitextra(ie)
            X = [X(1,:) - step(1); X];
            Y = [Y(1,:); Y];
            Z = [[result.extras.splitmeans(ie,:) result.extras.splitmeans(ie,end)]; Z];
        end
    end
    hold off;
    pcolor(X,Y,Z);
    if forcezero
        zlim = caxis;
        caxis([0 zlim(2)]);
    end
    colorbar;
    if(logx)
        set(gca,'Xscale','log','Yscale','log');
    end
end

  if showcounts && showplot && isfield(result,'y'); 
      ystep = mean([diff(result.y(1,:)) diff(result.y(:,1))'])/2;
      for ie = 1:size(result.x,2)
      for j = 1:length(result.x(:,ie,nc))
          if showpcolor
              text(result.x(j,ie)+step,result.y(j,ie)+ystep,sprintf('%d',result.n(j,ie)));
          elseif mksdf
              if showsdf(j,ie)
              text(times(end),result.sdfs{j,ie}(end),sprintf('%d',result.n(j,ie)),'color',colors{ie+addn});
              end
          else
              if logx
              text(result.x(j,ie)*stepratio,result.means(j,ie),sprintf('%d',result.n(j,ie)),'color',colors{ie+addn});
              else
              text(result.x(j,ie)+step,result.means(j,ie),sprintf('%d',result.n(j,ie)),'color',colors{ie+addn});
              end
          end
      end
  end
  end

  if ~isfield(result,'title')
      if exist('btype','var')
          result.title = MakeTitle(Expt,type,btype,stimtype);
      else
          result.title = MakeTitle(Expt,type,'00',stimtype);
      end
      if autonmin
          result.title = [result.title sprintf('nmin %0f',nmin)];
      end
      
      
      if showtitle & showplot
          title(result.title)
      end
  end


function val =SetLineval(bvals,xval)

if xval < -999 %% special interleave
   vmin = min(bvals(find(bvals > -1000)));
   step = mean(diff(bvals(find(bvals > -1000))));
end

if xval == -1003
    val = vmin - step/2;
elseif xval == -1004
    val = vmin - step;
elseif xval == -1001
    val = vmin - step/4;
elseif xval == -1002
    val = vmin - step/2;
else
    val = xval;
end

function str = MakeTitle(Expt,type,type2,stimtype)

if isfield(Expt.Header,'ExptFile')
    str = splitpath(Expt.Header.ExptFile);
else
    str = splitpath(Expt.Header.Name);
    if isfield(Expt.Header,'cellnumber')
        if isfield(Expt,'probes')
            str = strrep(str,'.mat',['.cell' num2str(Expt.Header.cellnumber) 'P' sprintf(':%d',unique([Expt.probes]))]);
        elseif Expt.Header.cellnumber ==0
            str = strrep(str,'.mat',[' MU P' sprintf(':%.1f',Expt.Header.probe)]);
        else
            str = strrep(str,'.mat',['.cell' num2str(Expt.Header.cellnumber) 'P' sprintf(':%.1f',Expt.Header.probe)]);
        end
    end
end
if ismember(stimtype,[7,17,19]) %% RDS/Sine
    str = sprintf('%s ic = %.2f',str,GetEval(Expt,'ic'));
elseif stimtype == 3
    str = [str sprintf(' Or = %.2f SF = %.2f TF = %.2f',GetEval(Expt,'or'),GetEval(Expt,'sf'),GetEval(Expt,'tf'))];
elseif stimtype == 18
    sf = GetEval(Expt,'sf');
    sfb = GetEval(Expt,'f2');
    str = [str sprintf(' SF = %.2f, +- %.2f ',sf,sfb-sf)];
elseif stimtype == 10
    str = [str sprintf(' SF = %.2f,%.2f TF = %.2f',GetEval(Expt,'sf'),GetEval(Expt,'f2'),GetEval(Expt,'tf'))];
elseif stimtype == 21 %image
%    str = [str sprintf('SF = %.2f',GetEval(Expt,'sf')) ' \pm ' sprintf('%.2f TF = %.2f  ',GetEval(Expt,'ob'),GetEval(Expt,'tf'))];
     str = [str sprintf(' SF = %.2f TF = %.2f  ',GetEval(Expt,'sf'),GetEval(Expt,'tf'))];
end
if strmatch(type2,'pi')
    or = GetEval(Expt,'or');
str = [str sprintf('%.1fx%.1f or%.0f',GetEval(Expt,'wi'),GetEval(Expt,'hi'),or)];
else
str = [str sprintf('%.1fx%.1f',GetEval(Expt,'wi'),GetEval(Expt,'hi'))];
end


function Fa = GetFa(Expt)

id = find([Expt.Trials.pi] > 0)
if isfield(Expt.Trials,'dfx')
    dfx = median([Expt.Trials(id).dfx] - [Expt.Trials(id).fx]);
else
    dfx = 0;
end
if isfield(Expt.Trials,'dfy')
    dfy = median([Expt.Trials(id).dfy] - [Expt.Trials(id).fy]);
else
    dfy = 0;
end
Fa = angle(dfx + i * dfy) .* 180/pi;
return


function str = val2str(val, type, stimtype,Expt,nfp)

if nargin < 5 | isempty(nfp)
    nfp = 2;
end

if(strmatch(type,'me'))
  if(val == -1)
    str = 'Left';
  elseif(val == 0)
    str = 'Binoc';
  else
    str = 'Right';
  end
elseif(strmatch(type,'xydir'))
   if val ==1
       str = 'Xpos';
   elseif val == 2
       str = 'Ypos';
   else
       str = sprintf('XYdir %.0f',val);
   end
elseif(strmatch(type,'sz'))
  str = sprintf('Size %.2f',val);
elseif(strmatch(type,'sd') & ismember(stimtype,[7 17 18]))
    if val == 0
        str = 'RDS';
    elseif val == 1
        str = 'R:grating L:rds';
    elseif val == -1
        str = 'R:rds L:grating';
    elseif val == 2
        str = 'Grating';
    else
        str = sprintf('%.2f',val);
    end
elseif(strmatch(type,'dq'))
    if(val == -1003)
        str = sprintf('sf%.2f',GetEval(Expt,'sf','mode'));
    elseif val == -1004
        str = sprintf('sf%.2f',GetEval(Expt,'f2','mode'));
    else
        str = sprintf('dq=%.2f',val);
    end
elseif(strmatch(type,'dp'))
    if(abs(val) == 1001)
        str = 'Right';
    elseif abs(val) == 1002
        str = 'Left';
    else
        str = sprintf('dp=%.2f',val);
    end
elseif(strmatch(type,'tf'))
  str = sprintf('tf=%.2f',val);
else
  str = sprintf('%s=%.*f',type,nfp,val);
end


    
function Expt = FillExpt(Expt,field)

if ~isfield(Expt.Trials,field)
  for j = 1:length(Expt.Trials)
      Expt.Trials(j).(field) = Expt.Stimvals.(field);
  end
else
  for j = 1:length(Expt.Trials)
      if isempty(Expt.Trials(j).(field))
      Expt.Trials(j).(field) = Expt.Stimvals.(field);
      end
  end
end

%
% sd -> CV (actually, vector length, is a Gaussian widh sd 28.64 degrees
function cv = sd2cv(sd)

cv = exp(-(sd.^2)./(2 .* 28.6462^2));

function sd = cv2sd(cv)

if cv > 1
    sd = 0;
else

sd = abs(sqrt(-log(cv) .* 2 .* 28.6462.^2));
end


function [xv, mindiff] = RemoveDuplicates(xvs, Ex)

xv = xvs(find(xvs > -1000));
step = mean(diff(sort(xv)));
mindiff = step/10;
if strmatch(Ex.Stimvals.et,{'dx' 'dy'}) %can be very small if dual increment
mindiff = range(xv)./200;
elseif strmatch(Ex.Stimvals.et,{'sO'}) %can be different for different thinor, 
mindiff = range(xv)./200;
else
mindiff = range(xv)./99;
end
if isempty(xv)
    return;
end

if isfield(Ex.Header,'setmindiff')
    mindiff = Ex.Header.setmindiff;
end
xv = sort(xv);
id = find(diff(xv) > mindiff);
if diff(size(id)) < 0
    id = id';
end
xv = xv([1 1+id]);



function [consistency, details] = CalcConsistency(Trials)
    ses = unique([Trials.se]);
    samedir = 0;
    diffdir = 0;
    srpts = [];
    for j = 1:length(ses)
        sid = find([Trials.se] == ses(j));
        srpts(j) = length(sid);
        for k = 2:length(sid) 
            if Trials(sid(k)).RespDir == Trials(sid(k-1)).RespDir
                samedir = samedir+1;
            else
                diffdir = diffdir+1;
            end
        end
    end
    consistency = samedir./(samedir+diffdir);
    details.seeds = ses;
    details.srpts = srpts;
    details.rptfraction = sum(srpts ==1)./length(ses);