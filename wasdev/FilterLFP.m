function [lfp, ft, removed]  = FilterLFP(lfp, samplerate,varargin)
%
%[lfp, ft] = FilterLFP(lfp, samplerate)
%Remove Mains freequencys from lfp. By default this is done in Fourier domain
%
%FilterLFP(lfp, samplerate,'cycavg') does this in the time domain,
%subtracting out the mean periodic variation at the mains period.
%


%mainsperiod = 166.6155;
mainsperiod = 166;
mainsf = 60;
%mainsperiod = 332;
id = 1;
j = 1;
navg = 0;
maxlen = 0;
setperiod = 166.624; %works well for first 50 trials of ruf648.rds.DT
setperiod = 166.72; %works well for ruf648.rds.DT with interp
ft = [];
removed = [];
while j < nargin -1
    if strncmpi(varargin{j},'period',4)
        j = j+1;
        if(varargin{j} > 0)
            setperiod = varargin{j};
        end
    elseif strncmpi(varargin{j},'length',3)
        j = j+1;
        if(varargin{j} > 0)
            maxlen = varargin{j};
        end
    elseif strncmpi(varargin{j},'freq',4)
        j = j+1;
        if(varargin{j} > 0)
            mainsf = varargin{j};
        end
    elseif strncmp(varargin{j},'FT',2)
        j = j+1;
        ft = varargin{j};
    elseif strncmpi(varargin{j},'cycavg',3)
        [lfp, ft, removed] = FilterLFPa(lfp, samplerate);
        return;
    elseif strncmpi(varargin{j},'avgb',3)
        [lfp, ft] = FilterLFPe(lfp, samplerate);
        return;
    end
    j = j+1;
end


%
%period = len/(10000 * samplerate)
%period of nth harmonic = len/(10000 * samplerate *n);
len = length(lfp);
aw=4;
nharm = 4;
if isempty(ft)
    ft = fft(lfp);
end
ft(1) = 0;
if(std(lfp) < 0.1)
    return;
end
for freq = [mainsf mainsf*2 mainsf*3]
    harmonic(1) = 1+round(len * freq/(10000 * samplerate));
    harmonic(2:(nharm+1)) = harmonic(1)+[1:nharm];
    harmonic = harmonic+1-floor(length(harmonic)/2); %+1 to igore DC, -2 to center 5 stim.
    idx = [(harmonic(1)-1-aw):(harmonic(1)-1) (harmonic(end)+1):(harmonic(end)+1+aw)];
    pwr(1) = mean(abs(ft(idx)));
    hpwr = abs(ft(harmonic));
    ft(harmonic) = ft(harmonic) .* pwr(1)./hpwr;
    rpts = length(ft)-(harmonic-2);
    ft(rpts) = ft(rpts) .* pwr(1)./hpwr;
end
lfp = real(ifft(ft)); %shoudl be real, but this deald with rounding.
%pwr(1) = mean(abs(ft(idx+1)));
%ft(harmonic+1) = mean(ft(idx+1));
%pwr(2) = mean(abs(ft(harmonic+1)));
%ft(harmonic+1) = ft(harmonic+1) .* pwr(1)/pwr(2);


function [lfp, ft, removed] = FilterLFPa(lfp, samplerate, start)

mainsperiod = 166.6155;
%mainsperiod = 166;
%mainsperiod = 332;
mainsperiod = 166;
id = 1;

%
% first caculate the mean periodic variation with period mainsperiod. 
% To avoid aliasing, build two averages with value added to bins either side
% of correct one. 
%
t = ((0:length(lfp)-1));
abins = 1+round(mod(t,floor(mainsperiod * samplerate)));
dbins = 1+ round(mod(t, ceil(mainsperiod * samplerate)));
for bin = 1:max(abins);
    amean(bin) = mean(lfp(find(abins == bin)));
end
for bin = 1:max(dbins);
    dmean(bin) = mean(lfp(find(dbins == bin)));
end
lfp = lfp - amean(abins)';
removed = amean;
%lfp = lfp - amean(abins)'/2;
%lfp = lfp - dmean(dbins)'/2;

ft = fft(lfp);

function [lfp, ft] = FilterLFPe(lfp, samplerate)

mainsperiod = 166.6155;
%mainsperiod = 166;
%mainsperiod = 332;
mainsperiod = 166.4;
id = 1;

%
% first caculate the mean periodic variation with period mainsperiod. 
% To avoid aliasing, build two averages with value added to bins either side
% of correct one. 
%
t = ((0:length(lfp)-1));
lfp = lfp - mean(lfp);
subsmp = 10;
abins = 1+round(mod(t*subsmp,floor(subsmp * mainsperiod * samplerate)));
for bin = 1:max(abins);
    nsmp(bin) = length(find(abins == bin));
    if(nsmp(bin) > 0)
        amean(bin) = mean(lfp(find(abins == bin)));
    end
end

lfp = lfp - amean(abins)';
ft = fft(lfp);

%ftfrq = (0:length(lfp)-1)/(length(lfp)/(10000 * samplerate));
