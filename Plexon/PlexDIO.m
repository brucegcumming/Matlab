function [adtimes, wtimes] = PlexDIO(loops, type, varargin)
%PlexDIO(nloops,type,...)  output pulses onto the DIO
% nloops in the number of pulses
% type 1 uses PL_DOOutputPulse line0, SetDuration line 0 - returns -10003)
% type 1 uses PL_DOOutputPulse - line1, SetDuration line 0 - returns -1
% type 3 uses PL_DOOutputPulse - line1, SetDuration line 1 both return -1
% otherwise uses PL_DOPulseBit - works

plclient = PL_InitClient(0);
pset = PL_GetPars(plclient);
[result, deviceNumbers, numBits, numLines] = PL_DOGetDigitalOutputInfo();
result =  PL_DOInitDevice(deviceNumbers(1), 1);
result = PL_DOClearAllBits(deviceNumbers(1));
bres = NaN; ares = NaN;
res = PL_DOSetLineMode(deviceNumbers(1), 1, 0)
res = PL_DOSetLineMode(deviceNumbers(1), 0, 0)

for j = 1:loops
    if type ==1
        ares = PL_DOSetPulseDuration(deviceNumbers(1), 0, 1000); % 1msec
        bres = PL_DOOutputPulse(deviceNumbers(1),0);
    elseif type ==2
        ares = PL_DOSetPulseDuration(deviceNumbers(1), 0, 1000); % 1msec
        bres = PL_DOOutputPulse(deviceNumbers(1), 1);
    elseif type ==3
        ares = PL_DOSetPulseDuration(deviceNumbers(1), 1, 1000); % 1msec
        bres = PL_DOOutputPulse(deviceNumbers(1),1);
    else
        ares =  PL_DOPulseBit(deviceNumbers(1), 1, 5); % 5msec
    end
    PL_DOSleep(200);
end
fprintf('Returns %d and %d\n',ares,bres);
PL_Close(plclient);
