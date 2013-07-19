function x = tryfit(Expt, varargin)

guess = [1 0.5 1 1 2 0 2];
%load('/bgc/bgc/data/rufus/924/ruf924.0.c1.grating.DPIC.mat')
%tryfit(Expt,'start',[1 0.5 1 1 4 pi 0])
% tryfit(Expt,'start',[0.5 0.5 0.8 1 6 pi -0.5]) gets close.
% Need to add - recifier. sqrt weighting, fminsearch
%
showguess = 0;
usesqrt = 1;
j = 1;
while j < nargin
    if strncmpi(varargin{j},'start',4)
        j = j+1;
        guess = varargin{j};
    elseif strncmpi(varargin{j},'showg',5)
        showguess = 1;
    end
    j = j+1;
end

res = PlotExpt(Expt);

Data.phases = res.x(:,1);
Data.contrasts(1) = max(res.y(1,:));
Data.contrasts(2) = min(abs(res.y(1,:)));
if showguess
    EvalFit(guess,Data,res,'plot');
end
options = optimset([]);

[x, rss] = fminsearch(@EvalFit, guess, options, Data, res);

EvalFit(x,Data,res,'phases',[-pi:0.1:pi],'plot');

function rss = EvalFit(A,Data,res, varargin)

%A(1) = ocularity -  1 = equal 0 = right eye 2 = left eye.
%A(2) = monoc divisor, low contrast
%A(3) = binoc divisor high-low combination
%A(4) = binoc divisior low both
%A(5) = Amplitude
%A(6) = Phase
%A(7) = mean

showplot = 0;
usesqrt = 0;
phases = Data.phases;
j = 1;
while j < nargin -2
    if strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'phases',5)
        j = j+1;
        phases = varargin{j};
    end
    j = j+1;
end


if A(1) < 1
    ocu(1) = A(1);
    ocu(2) = 1;
else
    ocu(1) = 1;
    ocu(2) = 1/A(1);
end
 
   amp(1) = ocu(1) .* Data.contrasts(1);
   amp(2) = ocu(2) .* Data.contrasts(2);
    [p,d,r,All] = grating(amp,[1 A(2) A(3)],'phases',phases,'dp',A(6));
    fitrates(:,1) = r;
 %
 % High contrast both eyes = no normalization
 %
   amp = ocu .* Data.contrasts(1);
    [p,d,r,All] = grating(amp,[1 1 1],'phases',phases,'dp',A(6));
    fitrates(:,4) = r;
  %low contrast both eyes  
    amp = ocu .* Data.contrasts(2);
     [p,d,r,All] = grating(amp,[A(2) A(2) A(4)],'phases',phases,'dp',A(6));
    fitrates(:,2) = r;
    
    amp(1) = ocu(1) .* Data.contrasts(2);
    amp(2) = ocu(2) .* Data.contrasts(1);
    [p,d,r,All] = grating(amp,[A(2) 1 A(3)],'phases',phases,'dp',A(6));
    fitrates(:,3) = r;
    
    fitrates = fitrates .* A(5) + A(7);
    fitrates(find(fitrates < 0)) = 0;
    if length(phases) == size(res.means,1)
        if usesqrt
            diffs = sqrt(fitrates)-sqrt(res.means);
        else
            diffs = fitrates-res.means;
        end
        rss = sum(sum(diffs.^2));
    end
    if showplot
        plot(phases,fitrates(:,1),'r');
        plot(phases,fitrates(:,2),'g');
        plot(phases,fitrates(:,3),'b');
        plot(phases,fitrates(:,4),'y');
    end






