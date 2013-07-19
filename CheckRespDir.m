function dir = CheckRespDir(Expt, type, varargin)
%returns the sign convention of respdir, based on the 
%monkey being right on average when there is a signal.
%positive answer means positive Resp matches positive signal
pid = find([Expt.Trials.(type)] > 0);
nid = find([Expt.Trials.(type)] < 0);
dir = mean([Expt.Trials(pid).RespDir])-mean([Expt.Trials(nid).RespDir]);