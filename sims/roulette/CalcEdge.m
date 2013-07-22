function [edge, redge, p] = CalcEdge(angles,varargin)


betw = 1;
resample = 0;
redge = [];
disttype = 0;
j = 1;
while j < nargin
    if strncmpi(varargin{j},'betw',3)
        j = j+1;
        betw = varargin{j};
    elseif strncmpi(varargin{j},'resample',3)
        j = j+1;
        resample = varargin{j};
    end
    j = j+1;
end
nrep = length(angles);
wid = find(abs(angles) < ((betw/2) * 2 * pi/37));
p = length(wid)/length(angles);
%     p * 36 - (1-p)  = p * 37 -1
% on average lose betw every go, gain 35 units each win.
% no, if count lose betw _every_ go, get 36 units each win
edge = (p * 36 - betw)/betw;

for j = 1:resample
    a = Bresample(angles);
    wid = find(abs(a) < ((betw/2) * 2 * pi/37));
    p = length(wid)/length(a);
%     p * 36 - (1-p)  = p * 37 -1
    redge(j) = (p * 36 - (1-p) * betw)/betw;
    %     p * 36 - (1-p)  = p * 37 -1
end
    