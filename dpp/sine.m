function resp = sine(params,x)
%Parameters are baseline, amplitude, phase.

resp =  params(1) +abs(params(2)) .* sin(x + params(3));
resp = [resp  zeros(length(resp),1)];
resp = max(resp,[],2);
