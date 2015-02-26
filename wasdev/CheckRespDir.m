function dir = CheckRespDir(Expt, type, varargin)
%returns the sign convention of respdir, based on the 
%monkey being right on average when there is a signal.
%positive answer means positive Resp matches positive signal
j = 1;
verbose = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'quiet',5)
        verbose = 0;
    elseif strncmpi(varargin{j},'verbose',5)
        verbose = 1;
    end
    j = j+1;
end
z = mean([Expt.Trials.(type)]);
pid = find([Expt.Trials.(type)] > z);
nid = find([Expt.Trials.(type)] < z);
presp = mean([Expt.Trials(pid).RespDir]);
sig(1) = mean([Expt.Trials(pid).(type)]);
sig(2) = mean([Expt.Trials(nid).(type)]);
nresp = mean([Expt.Trials(nid).RespDir]);
dir = presp-nresp;
if verbose
   fprintf('%s %.2f:%.3f Resp   %.2f:%.3f Resp %.2f\n',type, sig(1), presp, sig(2), nresp);
end