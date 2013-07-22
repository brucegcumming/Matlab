function ofrsum(varargin)

resptype = 'gauss';
sd = 4;
x = 1:10;

j = 1;
while j <= length(varargin)
    if strmatch(varargin{j},'parabolic')
    elseif strmatch(varargin{j},'sd')
        j = j+1;
        sd = varargin{j};
    end
    j = j+1;
end

if strmatch(resptype,'gauss')
    ri = exp(-(x-1).^2./sd.^2);
elseif strmatch(resptype,'parab')
    ri = 1 - 0.1*x.^2;
else
ri = 11-x;
end


for j = 1:1
    k = j*length(x);
    R(:,j) = cumsum(ri .* k)./x;
end
R(:,j+1) = cumsum(ri);
plot(x,R);
hold on;
plot(x,mean(R'),'k--');
plot(x,ri,'k-','linewidth',2);