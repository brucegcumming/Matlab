function res = CheckFullV(FullV, varargin)
% CheckFullV(FullV, varargin)
%Try triggering different ways. See if it is posible to figure out need for
%+ or - trigger

ispk = 22;
scales = 1;
spkrate = 100;
showth = 0;
filtershape = 0;
maxsize = NaN;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'maxsize',5)
        j = j+1;
        maxsize = varargin{j};
    elseif strncmpi(varargin{j},'odd',3)
        filtershape = 1;
    elseif strncmpi(varargin{j},'scale',4)
        j = j+1;
        scales = varargin{j};;
    elseif strncmpi(varargin{j},'tchan',4)
        j = j+1;
        ispk = varargin{j};
    end
    j = j+1;
end


if iscellstr(FullV)
    names = FullV;
    for j = 1:length(names)
        go = 1;
        if maxsize > 0
            d = dir(names{j});
            if d.bytes > maxsize
                fprintf('Ignoring %s size is %d\n',names{j},d.bytes);
                go = 0;
            end
        end
        if go
        fprintf('Loading %s\n',names{j});
        FullV = LoadFullV(names{j});
        res{j} = CheckFullV(FullV,varargin{:});
        res{j}.name = names{j};
        clear FullV; 
        end
    end
    return;
end
hold off;
means = [];
colors = mycolors;
duration = size(FullV.V,2).*FullV.samper;
nevents = duration .* spkrate;
for k = 1:length(scales)
    scale = abs(scales(k));
    if filtershape == 1 || scales(k) < 0 
        T = Gabor([0.1 0.35 * scale -pi/2 1 0 0],'npts',41);
    else
        a = Gauss(2 * scale,-20:20);
        b = Gauss([5*scale 8*scale],-20:20);
        T = a./sum(a)- (b./sum(b));
    end
    res.filters(k,:) = T;
T = T - mean(T);
T = T./std(T);
res.filters(k,:) = T;
for j = 1:length(ispk)
x(j,:) = conv(FullV.V(ispk(j),:),T);
sgn(j,:) = diff(sign(diff(x(j,:),1,2)),1,2);
pid = find(sgn(j,:) < 0)+1;
nid = find(sgn(j,:) > 0)+1;
px = x(j,pid);
t = prctile(px,100 - (100 .* nevents./length(px)));

nx = -x(j,nid);
tb = prctile(nx,100 - (100 .* nevents./length(nx)));
tlim(1) = min([px nx]);
tlim(2) = max([px nx]);
bins = linspace(tlim(1),tlim(2),100);
a = hist(px,bins)./length(px);
b = hist(nx,bins)./length(nx);
if length(scales) > length(ispk)
    c = colors{k};
else
    c = colors{j};
end
h = plot(bins,cumtrapz(bins,a)-cumtrapz(bins,b),'color',c);
if length(scales) < length(ispk) && k == 2
    set(h,'linestyle','--');
end
hold on;
if showth
plot([t t],get(gca,'ylim'));
plot([tb tb],get(gca,'ylim'),'r');
end
%plot(bins,a,'color',colors{k});
%plot(T,'color',colors{k});
a = cumsum(px)./length(px);
b = cumsum(nx)./length(nx);
hold on;
means(j,k) =  mean(px(px > t))-mean(nx(nx > tb));
end
end
if length(scales) >  length(ispk)
legend(num2str(scales'));
else
    legend(num2str(ispk'));
end

res.means = means;