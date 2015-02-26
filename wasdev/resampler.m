function [Kernel] = resampler(a,varargin)

%a is a vector of .mat psychfiles for ORBWCRC expts.
%resampler combines them into one expt. and does resampling.

rweight=0;
while j<=length(varargin)
    if strncmpi(j,'rweight',2)
        rweight=1;
    else
    end
    j=j+1;
end

Expt=CombineExpts({a.Kernel});
Expt = rmfield(Expt,'kernel')
if rweight==1
    [Kh Kernel] = BuildBWPsychRCa(Expt,'resample',1000,'rweight');
else
    [Kh Kernel] = BuildBWPsychRCa(Expt,'resample',1000);
end
Kernel.kernel=Kh;