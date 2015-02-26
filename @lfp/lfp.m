classdef lfp
%Builing/using LFP Files
%For Spike2 Data Files, uske MakeMat.s2s to generat raw lfp.mat files
%  Then LoadLFP(Expt) will add LFP data from the raw files into the Expt
% For Utah Recordings FullV2LFP.m coverts the FullV datafiles into LFP
%  Which can then be used by LoadLFP
    properties
        CurrentVersion = '1.0';
%1.14 Makes Sure that Spikes Values Matrix dims are same meaning, for all event counts        
    end
    methods (Static)        

    end
end