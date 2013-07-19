function rf = xytrf(varargin)

ori = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'ori',3)
        j = j+1;
        ori = varargin{j};
    end
end
x = -127:128;
y = -127:128;
t = -127:128; %ms

XY = Gabor([1 0.5 pi/2 1 0 0 ori 2]);
t = Gabor([0.1 0.5 pi/2 1 0 0]); %Odd symmetric Gabor for T

for j = length(t):-1:1
    rf(:,:,j) = XY .* t(j);
end