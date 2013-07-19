function OffFlash(varargin)
%OffFlash  Simulate LN responses to generate "Missing Flash" responses.


x = -50:150;
tk = Gauss([40 30 0.012],x)- Gauss(10,x);

npulse = 50;
pixedpulse = 0;

colors = 'brgkcm';
nr = 1;
nodc = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nodc',4)
        nodc = 1;
    end
    j = j+1;
end

hold off;
for period = [10 20];
stim = [zeros(1,100) round(mod(1:(npulse * period),period)./period) zeros(1,500)] .* -1;
if nodc
    stim = [zeros(1,100) round(mod(1:(npulse * period),period)./period)-0.5 zeros(1,500)] .* -1;
elseif fixedpulse
    stim = [zeros(1,100) round(mod(1:(npulse * period),period)./period) zeros(1,500)] .* -1;
else
    stim = [zeros(1,100) round(mod(1:(npulse * period),period)./period) zeros(1,500)] .* -1;
end
lresp = conv(tk,stim);
fid = find(stim < 0);
offtime = fid(end);
plot([1:length(lresp)]-offtime,lresp,'color',colors(nr));
hold on;
nr = nr+1;
end
