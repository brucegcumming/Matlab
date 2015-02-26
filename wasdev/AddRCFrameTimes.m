function Expt = AddRCFrameTimes(Expt, vars, varargin)
%Adds Start/End Times to an RC Expt
%Expt = AddRCFrameTimes(Expt, vars, varargin)
%vars is a cell array of field names that are defined for each frame
% This allows adjustments for dropped frames

framerpt = 1;

if isfield(Expt.Stimvals,'Fr') && Expt.Stimvals.Fr > 1
    framerpt = Expt.Stimvals.Fr;
end
frameperiod = Expt.Header.frameperiod;
timeshift = 15000;

for t = 1:length(Expt.Trials)
    T = Expt.Trials(t);
    if isfield(T,'Fr')
        frpt = T.Fr;
    else
        frpt = framerpt;
    end
    ev = T.(vars{1});
    if length(T.Start) == 1
        T.Start = timeshift + T.Start + [0:length(ev)-1]' .* frameperiod * frpt;
        T.End = T.Start + frameperiod * frpt;
        Expt.Trials(t) = T;
    end
end