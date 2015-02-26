%=================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GaborFit = FitSqrtGabor(x,y,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================================================================================================


% Fits a Gabor to the function y(x).
% NB the argument y is a cell array containing the different values obtained at each x.
% If a third argument is supplied, it is assumed to be an array containing the uncorrelated response
% If a fourth argument is supplied, it is assumed to be an initial-guess Gabor fit
% Takes sqrt of y and fits sqrt of a Gabor to it.
% so NB taing the the fit and squaring is not an exact fit to the raw means


if nargin==3
    sqrtyuncorr = sqrt(varargin{1});
    if ~isempty(sqrtyuncorr)
        uncorr = mean(varargin{1});
    end
elseif nargin==4
    sqrtyuncorr = sqrt(varargin{1});
    if ~isempty(sqrtyuncorr)
        uncorr = mean(varargin{1});
    end
    g = varargin{2};
    suppliedguessparams = [g.amp g.SD g.phase g.SF g.baseline g.offset];
end


% Work out mean
nx = length(x);
for j=1:nx
    meany(j) = mean(y{j});
    nrep(j) = length(y{j});
end


%%%%%%%%%%%
% Constraints on parameters - eg Nyquist limit, don't allow excessive amplitudes.
maxamp = range(meany)*2;
maxfreq = 0.25/mean(diff(x));
minoffset = min(x);
maxoffset = max(x);


% Convert cell array y into one long array
yy=[];
xx=[];
for j=1:nx
    xx = [xx x(j)*ones(size(y{j}))];
    yy = [yy y{j}];
end
sqrtyy = sqrt(yy);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get initial guesses


% i) he estimates basewline from end points
base = 0.5 * (meany(1) + meany(end));
if exist('uncorr')
    % If user gave us uncorr, may as well use this in inital estimate of baseline
    base = (meany(1) + meany(end) + uncorr)/3;
end


% ii) he estimates offset from centre of mass of squared diff from baseline:
A=0;B=0;
for j=1:nx
    A = A + sum((y{j}-base).^2) * x(j);
    B = B + sum((y{j}-base).^2) ;
end
offset = A/B;


% Do not allow offset to fall outside data range
if offset > max(x)
    offset = max(x);
elseif offset < min(x)
    offset = min(x)
end


% iii) Bruce gets freq from FT
DC = trapz(x,meany)/range(x);
freq = [0:1/256:1]*maxfreq;
firstmx=-Inf;
gotmax=0;
for jf=1:length(freq)
    FT(jf) = abs(trapz(x,(meany-DC) .* exp(2*pi*i*freq(jf).*x))) ;


    if FT(jf)>firstmx & gotmax==0
        firstmx = FT(jf);
        firstmxindx = jf;
    end
    if firstmx > 0 & FT(jf) < firstmx
        gotmax=1;
    end


end


exitflag=-1;
bestfval=Inf;
options = optimset('maxfunevals',1e7,'display','none','maxiter',1e3);
[mx,indx] = max(FT);
% If the user supplied an initial guess, use that; else work out a good starting point,
% first using the global maximum of the FT, then using the first peak of the FT
if indx == firstmxindx
    trials = 1;% First peak is also global maximum
else
    trials = [1 2]; % First local peak is not the same as global maximum, so try both in turn.
end
if exist('suppliedguessparams')
    trials = [trials 3];
end
randattempt = 0;
while exitflag<=0 & randattempt<100;
    randattempt = randattempt + 1;
    for trial=1:trials
        if trial==3
            guessparams = suppliedguessparams;
        else


            % Bruce's estimate of disparity freq is the position of the first local maximum in the FT of the DTC,
            % or the global maximum if this is 50% larger
            [mx,indx] = max(FT);
            if trial==1 % first try global max
                peakfreq = freq(indx);
            elseif trial==2 % then try first max
                peakfreq = freq(firstmxindx);
            end


            indx = min(find(FT<exp(-0.5)*mx & freq>peakfreq));
            freq_hisd = freq(indx);


            % iv) and SD ditto
            % Find where FT first drops belwo exp(-0.5) of peak
            indx = min(find(FT<exp(-0.5)*mx & freq>peakfreq));
            freq_hisd = freq(indx);
            if ~isempty(freq_hisd)
                SD = 1 / (2*pi*(freq_hisd-peakfreq));
            else
                SD = 1 / (2*pi*(max(freq)-peakfreq/2));
            end
            % v) he estimates phase and amplitude from Fourier spectrum:
            % y-base should be A exp(-u^2/2/SD^2) cos(2pi SF u + phase)
            % where u = x-offset.
            % FT of (y-base) at peak freq SF is approximately A*sqrt(pi/2)*SD * exp(-i*phase)
            % hence can get A.
            S = 0; C = 0;
            dx = x - offset;


            FTfit = trapz(dx, (meany-base) .* exp(2*pi*i.*peakfreq.*dx) );
            amp = abs(FTfit) / sqrt(pi/2) / SD;
            phase = -angle(FTfit);


            guessparams = [peakfreq, SD, phase, amp, offset, base];
        end
        % By this stage should have an initial guess, either from the user or worked out by the program
        % add noise to initial guess
        guessparams = guessparams + randn(size(guessparams)) .* 0.01 .* guessparams .* randattempt;
        % Don't allow negative amp, SD, peakfreq or base:
        for j=[1 2 4 5]; guessparams(j) = abs(guessparams(j)); end


        if exist('sqrtyuncorr') & ~isempty(sqrtyuncorr)
            [fitparams,fval,exitflag] = fminsearch(@SqrtGaborFitSSD2,guessparams,options,xx,sqrtyy,sqrtyuncorr,maxamp,maxfreq,minoffset,maxoffset);
        else
            [fitparams,fval,exitflag] = fminsearch(@SqrtGaborFitSSD2,guessparams,options,xx,sqrtyy,[],maxamp,maxfreq,minoffset,maxoffset);
        end
        if exitflag<=0 % try again
            if exist('sqrtyuncorr') & ~isempty(sqrtyuncorr)
                [fitparams,fval,exitflag] = fminsearch(@SqrtGaborFitSSD2,fitparams,options,xx,sqrtyy,sqrtyuncorr,maxamp,maxfreq,minoffset,maxoffset);
            else
                [fitparams,fval,exitflag] = fminsearch(@SqrtGaborFitSSD2,fitparams,options,xx,sqrtyy,[],maxamp,maxfreq,minoffset,maxoffset);
            end
        end


        if fval <= bestfval
            bestfitparams = fitparams;
            bestfval = fval;
        end
    end % next trial
    GaborFit.baseline = bestfitparams(6);
    GaborFit.amp = bestfitparams(4);
    GaborFit.SF = bestfitparams(1);
    GaborFit.phase = bestfitparams(3);
    GaborFit.SD = abs(bestfitparams(2));
    GaborFit.offset = bestfitparams(5);
    GaborFit.fval = fval;


end % next attempt


GaborFit.failed = 0;
if exitflag<=0 | bestfval==Inf
    disp('Fit failed!!!!')
    GaborFit.failed = 1;
end


