function Fit = FitSine(x, y, varargin)
%
% Fit = FitSine(x, y, varargin)
% Fits a halfwave rectified sinewave to y as a function of x.
% FitSine(params, x,'eval'), returns the fitted function
%
%Fit.answer = [baseline amplitude phase]
x = reshape(x,prod(size(x)),1);
y = reshape(y,prod(size(y)),1);

state.x = x;

dosqrt = 0;
halfrect = 0;
j = 1;
while j < nargin -1
    if strncmpi(varargin{j},'eval',4)
        Fit = sine(y,x);
        return;
    elseif strncmpi(varargin{j},'sqrt',4)
        dosqrt = 1;
    elseif strncmpi(varargin{j},'halfrect',4) % don't allow mean < 0
        state.halfrect = 1;
    end
    j = j+1;
end

x = x(find(~isnan(y)));
y = y(find(~isnan(y)));
state.x = x;
[guess(2), b] = famp(x,y,1/(2 * pi));
guess(1) = mean(y);
guess(3) = angle(b) + pi/2;
Fit.guess = guess;

if dosqrt
    Fit.answer = nlinfit(state,sqrt(y),@rootsine,guess);
else
    Fit.answer = nlinfit(state,y,@sine,guess);
end

Fit.dp = state.x;
Fit.sf = 1;
Fit.fitx = x;
Fit.answer(2) = abs(Fit.answer(2));
if isfield(state,'halfrect')
    Fit.answer(1) = abs(Fit.answer(1));
end
Fit.A = Fit.answer(2) .* cos(Fit.answer(3)-pi/2) + i * Fit.answer(2) .* sin(Fit.answer(3)-pi/2);
Fit.fitted = sine(Fit.answer,state);


function resp = sine(params,state)
%Parameters are baseline, amplitude, phase.

if isfield(state,'halfrect')
    params(1) = abs(params(1));
end
resp =  params(1) +abs(params(2)) .* sin(state.x + params(3));
resp = [resp  zeros(size(resp))];
resp = max(resp,[],2);


function resp = rootsine(params,x)
%Parameters are baseline, amplitude, phase.

resp =  params(1) +abs(params(2)) .* sin(x + params(3));
resp = [resp  zeros(size(resp))];
resp = sqrt(max(resp,[],2));

