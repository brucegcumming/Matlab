function RGC(varargin)
%RGC simulate simple effects of RGC RFs on stimuli, for Boris Expts
%

myparams = [];
simtype = 0;
heights = [10:20:500];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'oned',4)
        simtype =2;
    elseif strncmpi(varargin{j},'ofr',3)
        simtype = 3;
        j = j+1;
        ofrfit = varargin{j};
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            myparams = varargin{j};
        end
    elseif strncmpi(varargin{j},'alpha',3)
        j = j+1;
        alphas = varargin{j};
        ofrfit = repmat(ofrfit,length(alphas),1);
        ofrfit(:,1) = alphas;
    elseif strncmpi(varargin{j},'striph',6)
        simtype =1;
    end
    j = j+1;
end


if simtype == 1
    [x,y,exc]  = gauss2d(5,[-100:100]);
    [x,y,inh] = gauss2d(25,[-100:100]);
    exc = 2 * exc./sum(exc(:));
    inh = inh./sum(inh(:));
    rgc = exc-inh;
    period = 50;
    stimw = period*16;
    stimy = [1:stimw]-stimw/2;
    grating = sin(stimy.*(2 * pi)./ period);
    for j = 1:length(heights)
        id = round(mod(stimy,heights(j))./heights(j));
        stim(id ==0,:) = repmat(grating,sum(id==0),1);
        stim(id ==1,:) = repmat(-grating,sum(id==1),1);
        resp = conv2(stim,rgc,'valid');
        resps(j) = sum(resp(:).^2);
        mresps(j) = max(resp(:));
        xresps(:,j) = sum(resp.^2,2);
        id = length(resp)/2 + -heights(j)/2:(heights(j)/2)-1;
        aresps(j) = mean(xresps(id,j));
    end
    plot(heights./period,resps,'o-');
elseif simtype == 2
    heights = [20:20:200 400:200:5000];
    [x,y,exc]  = gauss2d(50,[-1000:1000]);
    [x,y,inh] = gauss2d(250,[-1000:1000]);
    exc = 2 * exc./sum(exc(:));
    inh = inh./sum(inh(:));
    rgc = exc-inh;
    period = 500;
    stimw = period*16;
    stimy = [1:size(rgc,1)]-size(rgc,1)/2
    grating = cos(stimy.*(2 * pi)./ period);
    prod = grating * rgc;
    fit = FitGauss(1:length(prod),prod);
    %for radii 5:1 and areas 2:1, this is a Gaussian with SD 0.1 cycles
    stimy = [1:stimw]-stimw/2;
    for j = 1:length(heights)
        stim = (2*round(mod(stimy,heights(j))./heights(j)))-1;
        resp = conv(stim,prod);
        resps(j) = max(resp(100:300));
        xresps(:,j) = resp;
        id = length(resp)/2 + [-heights(j)/2:(heights(j)/2)-1];
        aresps(j) = mean(resp(id).^2);
    end
    plot(heights./period,aresps,'o-');
elseif simtype == 3  %boris' equation.
   heights = [50:50:400 600:200:5000];
   period = 500;
   prod = Gauss(50, -1000:1000);
   stimw = period*16;
   stimy = [1:stimw]-stimw/2;
   for j = 1:length(heights)
       stim = (2*round(mod(stimy,heights(j))./heights(j)))-1;
       resp = conv(stim,prod);
       id = length(resp)/2 + [-heights(j)/2:(heights(j)/2)-1];
       M(j) = mean(resp(id).^2);
   end
   S = heights./period;
   for j = 1:size(ofrfit,1)
    x = ofrfit(j,:);
    ofr = (x(2).*M.^x(1))./(1 + x(3).*S);
    plot(S,ofr,'o-');
    hold on;
   end
    options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
    state.S = S;
    state.M = M;
    state.ofr = ofr;
    if length(myparams) == 3
        [ssd, myofr] = FitOFR(myparams, state);
        plot(S,myofr,'mo-');
        guess = myparams;
    else
        guess = x;
    end
    [fit,fval,exitflag, output] = fminsearch(@FitOFR,guess,options,state);
    [ssd, myofr] = FitOFR(fit, state);
    plot(S,myofr,'ro-');
else
    exc  = Gauss(5,[-200:200]);
    inh = Gauss(25,[-200:200]);
    exc = 2 * exc./sum(exc);
    inh = inh./sum(inh);
    rgc = exc-inh;
    ft = fft(rgc);
    plot(abs(ft));
[a,b] = max(abs(ft(1:200)));
end

function [SSD, myofr ] = FitOFR(x, state)

myofr = (x(2).*state.M.^x(1))./(1 + x(3).*state.M.^10)
diffs = myofr-state.ofr;
SSD = sum(diffs.^2);
