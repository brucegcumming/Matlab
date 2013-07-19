function [adtimes, wtimes] = ReadPlex(varargin)


plclient = PL_InitClient(0);
pset = PL_GetPars(plclient);
plotwd = 0;
plotad = 2;

go = 1;
loops = 1;
nloops = 100;

if isnumeric(varargin{1})
    nloops = varargin{1};
    j = 2;
else
    j = 1;
end

while j <= length(varargin)
    if strncmpi(varargin{j},'plotwave',6)
        plotwd = 1;
    elseif strncmpi(varargin{j},'plotad',6)
        plotad = 1;
    elseif strncmpi(varargin{j},'noplot',6)
        plotad = 0; plotwd = 0;
    end
    j = j+1;
end
wtimes = [];
alladtimes = [];
allad = [];

[result, deviceNumbers, numBits, numLines] = PL_DOGetDigitalOutputInfo();
result =  PL_DOInitDevice(deviceNumbers(1), 1);
result = PL_DOClearAllBits(deviceNumbers(1));


while(go)
    if rand(1,1) < 0.1
        a = PL_DOSetPulseDuration(deviceNumbers(1), 1, 10000);
        a = PL_DOOutputPulse(deviceNumbers(1),1);
%  a=  PL_DOPulseBit(deviceNumbers(1), 1, 5); % 5 msec
    end
    PL_WaitForServer(plclient,1000);
    [n,t,d] = PL_GetAD(plclient);
    [wn, wt, wd] = PL_GetWF(plclient);
    [a, b] = PL_GetTS(plclient);
    if a > 0
        size(b)
    end
    adtimes(loops) = t;
    at = t + [1:size(d,1)]./pset(8);
    alladtimes = [alladtimes at];
    if wn > 0 
        if plotwd
            plot(wd');
        end
        wtimes = [wtimes wt'];
    end
    if n > 0 
        allad = [allad; d];
        if plotad == 2
            plot(allad);
        elseif plotad == 1
            plot(d);
        end
    loops = loops+1;
    elseif wn > 0
        loops = loops+1;
    end
    if loops > nloops
        go = 0;
    end
%    plot(wd); 
    drawnow;
end
PL_Close(plclient);
diff(adtimes);