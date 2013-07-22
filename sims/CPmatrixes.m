function result = CPmatrixes(varargin)


x = 0:5:360;
mode = 2; %1 = corr = cos(diff), 2 = sign(diff) (shadlen), 3 = diff (linear)
task = 1;
poolw = 20;
poolo = 0;
for j = 1:length(varargin)
    if strncmpi(varargin{j},'discrim',5)
        poolw = 90;
        task = 1;
    elseif strncmpi(varargin{j},'fine',4)
        task = 2;
        poolw = 20;
    elseif strncmpi(varargin{j},'mkfig2',6)
        a = cpmatrixes('linear','poolw',1);
        b = cpmatrixes('linear','poolw',45);
        c = cpmatrixes('cosine','poolw',1);
        d = cpmatrixes('cosine','poolw',45);
        GetFigure('CPOri');
        subplot(1,1,1);
        hold off;
        plot(a.x,a.cdiff(:,1),'b--');
        hold on;
        plot(b.x,b.cdiff(:,1),'b');
        plot(c.x,c.cdiff(:,1),'r--');
        plot(d.x,d.cdiff(:,1),'r');
        figify(gcf,gca,'fontsize',12,'fontname','Ariel');
        set(gca,'xlim',[0 360],'ytick',[],'xtick',[0 90 180 270 360],'xticklabel',{'-180' '-90' '0' '90' '180'});
        return;
    elseif strncmpi(varargin{j},'makefigures',7)
        GetFigure('Shadlen10')
        CPmatrixes('shadlen');
        GetFigure('Shadlen45');
        CPmatrixes('shadlen','poolwidth',45);
        
        GetFigure('linear45');
        CPmatrixes('linear','poolwidth',45);
        GetFigure('linear10');
        CPmatrixes('linear','poolwidth',10);
        GetFigure('cos45');
        CPmatrixes('cosine','poolwidth',45);
        GetFigure('cos10');
        CPmatrixes('cosine','poolwidth',10);
        return;
    elseif strncmpi(varargin{j},'shadlen',5)
        mode = 2;
    elseif strncmpi(varargin{j},'linear',6)
        mode = 3;
    elseif strncmpi(varargin{j},'cosine',5)
        mode = 1;
    elseif strncmpi(varargin{j},'pooloffset',5)
        j = j+1;
        poolo = varargin{j};
    elseif strncmpi(varargin{j},'poolwidth',5)
        j = j+1;
        poolw = varargin{j};
    end
    j = j+1;
end

if task == 2
    poola  = find(x >90+poolo & x <= 90++poolo+poolw);
    poolb = find(x >= 90-poolw-poolo & x < 90-poolo);
else
poola  = find(x > 90 - poolw & x < 90 + poolw);
poolb = find(x > 270-poolw & x < 270+poolw);
end

cmax = 0.1;
for j = 1:length(x)
    d = x(j) -x;
    if mode == 1
        C(:,j) = [1 + sin(d * pi/180)] .* cmax/2;
    elseif mode == 2
        C(:,j) = [1 + sign(sin(d * pi/180))] .* cmax/2;
    elseif mode == 3
        C(:,j) =  [0.5+ asin(sin(d * pi/180))./pi] .*cmax;
    end
end
subplot(2,1,1);
hold off;
imagesc(x,x,C);
set(gca,'ydir','normal');
hold on;
a = min(x(poola));
b = max(x(poola));
plot([a b b  a a],[a a b b a],'w');
a = min(x(poolb));
b = max(x(poolb));
plot([a b b  a a],[a a b b a],'w');
for j = 1:size(C,1)
    corr(j,1) = mean(C(poola,j));
    corr(j,2) = mean(C(poolb,j));
end
subplot(2,1,2);
hold off;
plot(x,corr);
hold on;
plot(x,corr(:,1)-corr(:,2),'r');
result.x = x;
result.corr = corr;
result.CMat = C;
result.cdiff = corr(:,1) - corr(:,2);