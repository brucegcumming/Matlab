classdef expt
    properties
        Version = '1.1'
    end
    
    
    
    methods (Static)
        Expt = InsertTrials(Expt, Trials, varargin)
        Expt = CountEvents(Expt, t, latency, varargin)
        Expt = AddComments(Expt, varargin)
        [type, details] = Type(name, varargin);
        X = TrialVals(Expt,  varargin)
        spk = spktimes(E, varargin)
        [xc, details] = xcorr(E, varargin)
        E = CheckFields(E,varargin)
        Expt = fix(Expt, type, varargin)
        rs = Resample(Expt, nresample, varargin);
        cp(src, tgt, varargin);
        [Expt, varargout] = SelectTrials(Expt, varargin)
end
end
