function Expts = SetTrialOffsets(Expts, varargin)
%SetTrialOffsets(Expts, varargin) ensure Trials.Trial increases
%monotonically.  This
% way .id can be left alone as it matches the serial output
% Trials.Trial is iimportant for making .Trial consecutive, and used for
% Header.BlockStart

Trialids = [];
for nexp = 1:length(Expts)
    E = Expts{nexp};
    if ~isempty(E)
    if isempty(Trialids)
        trialoffset = 0;
    else
        trialoffset = max(Trialids) - E.Trials(1).Trial+1;
        if trialoffset < 0
            trialoffset = 0;
        end
    end
    if isfield(E.Header,'trialoffset') && E.Header.trialoffset ~= trialoffset && E.Header.trialoffset < 0
        cprintf('red','Changing trial offset to %d (was %d) for %s\n',trialoffset,E.Header.trialoffset,GetEval(E,'name'));
    end
    E.Header.trialoffset = trialoffset;
    for j = 1:length(Expts{nexp}.Trials)
        E.Trials(j).Trial = E.Trials(j).Trial+trialoffset;
    end
    newt = [E.Trials.Trial];
    Trialids = [Trialids newt];
    Expts{nexp} = E;
    end
    trialoffsets(nexp) = trialoffset;
end
find(diff(Trialids) < 0);