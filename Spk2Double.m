function S = Spk2Doulbe(S)
%Spk2Doulbe(S)Converts ints to double in Spks structs
if isfield(S,'maxv') && isinteger(S.values)
    S.values = double(S.values) .* S.maxv./S.maxint;
end
if isfield(S,'xvalues') && isinteger(S.xvalues)
    S.xvalues = double(S.xvalues) .* S.xmaxv./S.maxint;
end
