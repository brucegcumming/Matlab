function Trials = MakeDummyTrials(varargin)
%Trials = MakeDummyTrials(...)
%
%


emtype = 2;
j = 1;
while j < nargin -1
    if strncmpi(varargin{j},'emtype',2)
        j = j+1;
        emtype = varargin{j};
    end
    j = j+1;
end

for j = 1:10;
    if emtype == 1
        x = rand(1600,1) - 0.5;
        Trials(j).Eyevals.lv = cumsum(x);
    else
        x = -800:800;
        Trials(j).Eyevals.lv = [x > 0]' - 0.5;
    end
end