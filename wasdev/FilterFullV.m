function FullV = FilterFullV(FullV,varargin)


boxw = 40;
sd = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'box',3)
        j = j+1;
        boxw = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sd = varargin{j};
    end
    j = j+1;
end
if sd > 0
   w = round(sd);
   kernel = Gauss(sd, -w*3:w*3);
   kernel = kernel./sum(kernel);
   o = floor(length(kernel)/2);
else
    kernel = ones(1,boxw)./boxw;
    o = floor(boxw/2);
end
sm = conv(FullV.V,kernel);
sm = sm(o:o+length(FullV.V)-1);
FullV.V = FullV.V -sm;