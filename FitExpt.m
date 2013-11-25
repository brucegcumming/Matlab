function fit = FitExpt(dat, varargin)
%FitExpt(Expt....)
%vanilla fits for RF measures
%Expt is a result output from PlotExpt, not a standard expt structure. 

fit = [];
showplot = 0;
fitargs = {};
fitted = [];
vonmises = 0;
specialfit = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plotfit',6)
        showplot = 2;
    elseif strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'vonmises',4)
        vonmises = 1;
    else
        fitargs = {fitargs{:} varargin{j}};
    end
     j = j+1;
end

if isempty(dat) 
    return;
end

if ~isfield(dat,'x')
    cprintf('red','No X in %s\n',dat.name);
    return;
end
    if ~isfield(dat,'handles')
        dat.handles = [];
    end

if isfield(dat,'type') 
    for j = 1:length(dat.type)
        if isempty(dat.type{j})
            dat.type{j} = 'e0';
        end
    end
end

if isfield(dat,'Data')
T = dat.Data.Trials;
else 
    T = [];
end

if isfield(dat,'bestdelay') %its an RC expt
    if isnan(dat.bestdelay)
        return;
    elseif strmatch(dat.type{1},{'Op' 'Pp'})
        x = dat.x(:,1);
        y = sum(dat.y(:,1,dat.bestdelay),2);
        if strmatch(dat.type{2},'ph')
            black = dat.y(:,1,dat.bestdelay);
            white = dat.y(:,2,dat.bestdelay);
            bnp = dat.sdfs.n(:,1);
            wnp = dat.sdfs.n(:,2);
        else
        black = [];
        white = [];
        bnp = [];
        wnp = [];
        end
        id = find(dat.sdfs.extraval == -1009);
        np = sum(dat.sdfs.n,2);
        bv = mean(dat.vars(1:10));
        sv = std(dat.vars(1:10));
        if id
           [a, t] = min(abs(dat.times./10-dat.delaysamples(dat.bestdelay)));
           x = [x; Inf];
           y = [y; dat.sdfs.extras{id}.sdf(t)];
           black = [black; dat.sdfs.extras{id}.sdf(t)];
           white = [white; dat.sdfs.extras{id}.sdf(t)];
           np = [np; dat.sdfs.extras{id}.n];
           bnp = [bnp; dat.sdfs.extras{id}.n];
           wnp = [wnp; dat.sdfs.extras{id}.n];
        end
        fit = FitGauss(x,y,'nreps',np,'freebase', fitargs{:});
        fit.resp = y;
        fit.x = x; %anybody need this too?  Tuning curve fits use x...
        fit.rcvar = (max(dat.vars)-bv)/sv;
        if length(black) > 1
            fit.bfit = FitGauss(x,black,'nreps',bnp,'freebase', fitargs{:});
            fit.bfit.respsd = std(black);
            fit.bfit.max = max(black);
            fit.wfit = FitGauss(x,white,'nreps',wnp,'freebase', fitargs{:});
            fit.wfit.respsd = std(white);
            fit.wfit.max = max(white);
             bvar = var(cell2mat(dat.sdfs.s(:,1)'),0,2);
             wvar = var(cell2mat(dat.sdfs.s(:,2)'),0,2);
            fit.btfit = FitGauss(dat.times,bvar','freebase', fitargs{:});
            fit.wtfit = FitGauss(dat.times,wvar','freebase', fitargs{:});
            fit.toverlap = (fit.btfit.mean-fit.wtfit.mean)./sqrt(mean([fit.btfit.sd^2 fit.wtfit.sd^2]));
        end
        id = find(~isnan(fit.x) & ~isinf(fit.x));
        fit.xv = linspace(min(fit.x(id)),max(fit.x(id)));
        fit.fitcurve = FitGauss(fit.xv,fit.params,fitargs{:},'eval');
        f = dat.type{1};
        fit = CopyFields(fit, dat.Data.Stimvals,{'rOp' 'rPp' 'Op' 'Pp'});
        if isfield(T,'xo') || isfield(T,'yo')
            xo = [];
            yo = [];
            for j = 1:length(T)
                if isfield(T,'xo')
                    xo(j) = T(j).xo;
                end
                if isfield(T,'yo')
                    yo(j) = T(j).yo;
                end
                op(j) = T(j).(f)(end);
            end
            id = find(op > -1000);
            fit.(['stim' f]) = mean(op(id));
            if isempty(xo)
                fit.opdir(1) = 0;
                fit.xo = dat.Data.Stimvals.xo;
            else
                fit.xo = mean(xo(id));
            xf = polyfit(op(id),xo(id),1);
            fit.opdir(1) = xf(1);
            end
            if isempty(yo)
                fit.opdir(2) = 0;
                fit.yo = dat.Data.Stimvals.yo;
            else
                fit.yo = mean(yo(id));
            xf = polyfit(op(id),yo(id),1);
            fit.opdir(2) = xf(1);
            end
        end
            
    elseif strmatch(dat.type{1},{'or'})
        x = dat.x(:,1);
        y = sum(dat.y(:,:,dat.bestdelay),2);
        np = sum(dat.sdfs.n,2);
        period = 180;
        bl = find(dat.sdfs.extraval == -1009);
        if bl && dat.sdfs.extras{bl}.n > 0
            x = [x; Inf];
            y = [y; dat.sdfs.extras{bl}.sdf(dat.delaysamples(dat.bestdelay))];
            np = [np; dat.sdfs.extras{bl}.n];
        end
        
            if vonmises
                fitargs = {fitargs{:} 'posbase'};
                fit = FitGauss(x,y,'nreps',np,fitargs{:},'periodwrap',period,'posamp');
                fit.x = x;
                id = find(~isnan(fit.x) & ~isinf(fit.x));
                fit.xv = linspace(min(fit.x(id)),max(fit.x(id)));
                fit.fitcurve = FitGauss(fit.xv,fit.params,fitargs{:},'periodwrap',period,'eval');
            else
                fitargs = {fitargs{:} 'freebase'};
                fit = FitGauss(x,y,'nreps',np,fitargs{:},'period',period);
                fit.x = x;
                id = find(~isnan(fit.x));
                fit.xv = linspace(min(fit.x(id)),max(fit.x(id)));
                fit.fitcurve = FitGauss(fit.xv,fit.params,fitargs{:},'period',period,'eval');
            end
            fit.peak = PeakFit(fit);
    end
elseif strmatch(dat.type{1},{'Op' 'Pp'}) %not RC expt
    x = reshape(dat.x,1,prod(size(dat.x)));
    y = reshape(dat.means,1,prod(size(dat.means)));
    np = reshape(dat.n,1,prod(size(dat.n)));
    bl = strmatch('Blank',dat.extras.label);
    if bl
        x = [x Inf];
        y = [y dat.extras.means(bl)];
        np = [np dat.extras.n(bl)];
        fitargs = {fitargs{:} 'freebase'};
    end
    fit = FitGauss(x,y,'nreps',np,fitargs{:});
    fit.x = x;
    fit.resp = y;
    fit.ve = 1 - fit.rss./(var(y) .* length(y)-1);
    theta = pi/2 - (GetEval(dat.Data,'Ro') * pi/180);
    y = GetEval(dat.Data,'Pp');
    xx = x .* cos(theta) + y * sin(theta);
    yy = y * cos(theta) - x * sin(theta);
    id = find(~isnan(fit.x) & ~isinf(fit.x));
    fit.xv = linspace(min(fit.x(id)),max(fit.x(id)));
    fit.fitcurve = FitGauss(fit.xv,fit.params,fitargs{:},'eval');
    fit.peak = PeakFit(fit);
    f = dat.type{1};
    id = find([dat.Data.Trials.(f)] > -1000);
    if isfield(dat.Data.Trials,'xo')
        xf = polyfit([dat.Data.Trials(id).(f)],[dat.Data.Trials(id).xo],1);
        fit.opdir(1) = xf(1);
    end
    if isfield(dat.Data.Trials,'yo')
        yo = [dat.Data.Trials(id).yo];
        if length(yo) == length(id)
            xf = polyfit([dat.Data.Trials(id).(f)],yo,1);
            fit.opdir(2) = xf(1);
        end
    end
        

elseif strmatch(dat.type{1},{'stimxy'}) %not RC expt
    amps = squeeze(var(dat.means));
    [a,b] = max(amps(:));
    [xi,yi] = ind2sub(size(amps),b);
    ors = unique(dat.z);
    thinor = squeeze(dat.thinor(1,:,:));
    bestor = dat.thinor(1, xi,yi);
    [a,b] = min(abs(cos((bestor - thinor(:)) .* pi/180)));
    [ai,bi] = ind2sub(size(amps),b);
    orthor = dat.thinor(1,ai,bi);
    id = find(thinor == orthor);
    [a,b] = max(amps(id));
    [ai,bi] = ind2sub(size(amps),id(b));
    x = dat.x(:,xi, yi)';
    y = dat.means(:,xi, yi)';
    np = dat.n(:,xi, yi)';
    z = dat.x(:,ai, bi)';
    zy = dat.means(:,ai, bi)';
    znp = dat.n(:,ai, bi)';
    bl = strmatch('Blank',dat.extras.label);
    if bl
        x = [x Inf];
        y = [y dat.extras.means(bl)];
        np = [np dat.extras.n(bl)];
        z = [z Inf];
        zy = [zy dat.extras.means(bl)];
        znp = [znp dat.extras.n(bl)];
        fitargs = {fitargs{:} 'freebase'};
    end
    fit = FitGauss(x,y,'nreps',np,fitargs{:});   
    fit(2) = FitGauss(z,zy,'nreps',znp,fitargs{:});
    fit(1).x = x;
    fit(2).x = z;
    fit(1).resp = y;
    fit(2).resp = zy;
    id = find(~isnan(fit(1).x) & ~isinf(fit(1).x));
    fit(1).xv = linspace(min(fit(1).x(id)),max(fit(1).x(id)));
    fit(1).fitcurve = FitGauss(fit(1).xv,fit(1).params,fitargs{:},'eval');
    id = find(~isnan(fit(2).x) & ~isinf(fit(2).x));
    fit(2).xv = linspace(min(fit(2).x(id)),max(fit(2).x(id)));
    fit(2).fitcurve = FitGauss(fit(2).xv,fit(2).params,fitargs{:},'eval');
    sina = sin(bestor*pi/180);
    if sina > 0.9 %90 degrees
        fit(1).RFpos(1) = fit(1).mean;
        fit(1).RFpos(2) = fit(2).mean;
        fit(1).color = 'r';
        fit(2).color = 'b';
    elseif sina > 0.6 %45 degrees
        fit(1).RFpos(1) = (fit(1).mean+fit(2).mean)./sqrt(2);
        fit(1).RFpos(2) = (fit(2).mean-fit(1).mean)./sqrt(2);
        fit(1).color = 'm';
        fit(2).color = 'g';
    elseif sina > -0.1 %0 degrees
        fit(1).RFpos(1) = fit(2).mean;
        fit(1).RFpos(2) = fit(1).mean;
        fit(1).color = 'b';
        fit(2).color = 'r';
    elseif sina > -0.8 %-45 degrees
        fit(1).RFpos(1) = (fit(1).mean+fit(2).mean)./sqrt(2);
        fit(1).RFpos(2) = (fit(1).mean-fit(2).mean)./sqrt(2);
        fit(1).color = 'g';
        fit(2).color = 'm';
    else
        fit(1).RFpos(1) = fit(1).mean;
        fit(1).RFpos(2) = fit(2).mean;
    end
    fit(1).fitstr = sprintf('RF %.2f,%.2f',fit(1).RFpos(1),fit(1).RFpos(2));
        
    
elseif strmatch(dat.type{1},{'or'}) & length(dat.type) > 1 & strmatch(dat.type{2},{'me' 'ob'})  %%nned period 360 or 180 
    x = dat.x;
    y = dat.means;
    np = dat.n;
    fitargs = {fitargs{:} 'posbase'};
    
    for j = 1:size(dat.x,2)
        specialfit = 0;
        x = dat.x(:,j);
        y = dat.means(:,j);
        id = find(~isinf(x) & ~isnan(y));
        x = x(id);
        y = dat.means(id,j);
        np = dat.n(id,j);
        if range(x) > 270
            r(1) = MeanVector(y,x);
            r(2) = MeanVector(y,x,'double');
            if abs(r(2)) > abs(r(1)) %not directional but wide ragne
                period = 180;
                specialfit = 1;
            else
                period = 360;
            end
        else
            period = 180;
        end
        bl = strmatch('Blank',dat.extras.label);
        if bl
            x = [x; Inf];
            y = [y; dat.extras.means(bl)];
            np = [np; dat.extras.n(bl)];
            fitargs = {fitargs{:} 'posbase'};
        end
        id = find(~isnan(y));
        if specialfit
            fit = FitOriTuning(x(id),y(id),'nreps',np(id),fitargs{:});
            fit.type = 'FitOriTuning';
        elseif vonmises
            fit = FitGauss(x(id),y(id),'nreps',np(id),fitargs{:},'periodwrap',period);
            fit.type = 'FitGauss';
        else
            fit = FitGauss(x(id),y(id),'nreps',np(id),fitargs{:},'period',period);
            fit.type = 'FitGauss';
        end
        fit.x = x;
        fit.resp = y;
        id = find(~isnan(fit.x) & ~isinf(fit.x));
        fit.xv = linspace(min(fit.x(id)),max(fit.x(id)));
        if specialfit
            fit.fitcurve = FitOriTuning(fit.xv,fit.params,fitargs{:},'eval');
        else
            fit.fitcurve = FitGauss(fit.xv,fit.params,fitargs{:},fit,'eval');
        end
        fit.respvar = std(y);
        fit.peak = PeakFit(fit);
        fits(j) = fit;
    end
    fit = fits;
elseif strmatch(dat.type{1},{'or'})  %%nned period 360 or 180 
    x = reshape(dat.x,1,prod(size(dat.x)));
    y = reshape(dat.means,1,prod(size(dat.means)));
    np = reshape(dat.n,1,prod(size(dat.n)));
    
    if range(x) > 270
        period = 360;
    else
        period = 180;
    end
    bl = strmatch('Blank',dat.extras.label);
    if bl && dat.extras.n(bl) > 0
        x = [x Inf];
        y = [y dat.extras.means(bl)];
        np = [np dat.extras.n(bl)];
    end

    if vonmises
        fitargs = {fitargs{:} 'posbase'};
        fit = FitGauss(x,y,'nreps',np,fitargs{:},'periodwrap',period,'posamp');
        fit.x = x;
        id = find(~isnan(fit.x) & ~isinf(fit.x));
        fit.xv = linspace(min(fit.x(id)),max(fit.x(id)));
        fit.fitcurve = FitGauss(fit.xv,fit.params,fitargs{:},'periodwrap',period,'eval');
    else
        fitargs = {fitargs{:} 'posbase'};
        fit = FitGauss(x,y,'nreps',np,fitargs{:},'period',period);
        fit.x = x;
        id = find(~isnan(fit.x) & ~isinf(fit.x));
        fit.xv = linspace(min(fit.x(id)),max(fit.x(id)));
        fit.fitcurve = FitGauss(fit.xv,fit.params,fitargs{:},'period',period,'eval');
    end
        fit.peak = PeakFit(fit);
        fit.resp = y;
   
elseif strmatch(dat.type{1},{'dx'})
    for j = 1:length(dat.x)
        x(j) = dat.x(j);
        y{j} = dat.counts{j};
    end
    fit = FitGabor(x,y);
    fit.resp = y;
else
    fit = [];
end

if ~isempty(fit)
if ~isfield(fit,'xv')
    j = 1;
end
if isfield(dat.Data.Stimvals,'ve')
    [fit.ve] = deal(dat.Data.Stimvals.ve);
end

if isfield(dat.Data.Trials,'xo')
    xos = cat(1,dat.Data.Trials.xo);
    xo = mean(xos(xos > -1000));
    [fit.xo] = deal(xo);
end
if isfield(dat.Data.Trials,'yo')
    xos = cat(1,dat.Data.Trials.yo);
    yo = mean(xos(xos > -1000));
    [fit.yo] = deal(yo);
end
if showplot == 1  
        id = find(~isnan(x) & ~isinf(x));
        plot(fit.x(id),fit.fitted(id));
    elseif showplot == 2
        if isempty(fitted)
            for j = 1:length(fit)
                id = find(~isnan(x) & ~isinf(x));
                xv = linspace(min(fit(j).x(id)),max(fit(j).x(id)));
                h = plot(fit(j).xv,fit(j).fitcurve);
                if isfield(fit,'color')
                    set(h,'color',fit(j).color);
                elseif length(dat.handles) >= j && ishandle(dat.handles(j))
                    set(h,'color',get(dat.handles(j),'color'));
                end
            end
        else
            for j = 1:size(fitted,1)
            plot(xv, fitted(j,:));
            end
        end
    end
end


function peak = PeakFit(fit)
if isfield(fit.state,'period')
        if fit.amp < 0
            peak = fit.mean-fit.state.period/2;
            if peak < min(fit.x)
                peak = peak + fit.state.period;
            end
        else
            peak = fit.mean;
        end
else
    peak = fit.mean;
end
