function SkewResample(varargin)
%SkewResample 
%Simulate Bootstrapping the difference of two means drawn from the same
%skewed distribution, but with diffferent N
%Plots an image of skewness of resampled differences, for each N1,N2 pair
%Press on squares in the image to plot the histogram of resampled
%differnces. Shift/Cntrl hit for separate distributions, skewness histogram
%Skew generated with exp(randn()/skewfactor)), so larger skewfactor is less
%skewed
%SkewResample('skewfactor',k,...) sets skewfactor to k
%SkewResample('nruns',N,...) sets the number of loops through the
%sample-resample loop (setting total number of trials for distribution);

skewfactor = 1.5; %smaller numbers make more skew
nsmpls = [5 10 20 40 80 160]; 
nruns = 100; %number of times to repeat simulation for smoother results

j = 1;
while j < length(varargin)
    if strncmpi(varargin{j},'fine',4)
        nsmpls = [5 10 20 40 80 160];
        nruns = 1000; %number of times to repeat simulation for smoother results
    elseif strncmpi(varargin{j},'skewfactor',4)
        j = j+1;
        skewfactor = varargin{j};
    elseif strncmpi(varargin{j},'quick',4)
        nsmpls = [10  40 160];
        nruns = 20; %number of times to repeat simulation for smoother results
    elseif strncmpi(varargin{j},'skewfactor',4)
        j = j+1;
        skewfactor = varargin{j};
    end
    j = j+1;
end


allx = [];
for j = 1:length(nsmpls)
    for k = 1:j
        nsmpl = [nsmpls(j) nsmpls(k)];
        alla = [];
        allb = [];
        for n = 1:nruns
            x = exp(randn(sum(nsmpl),1)./skewfactor);
            allx = [allx x'];
            X = Bresample(x,1000);
            a = mean(X(:,1:nsmpl(1)),2);
            b = mean(X(:,nsmpl(1)+1:end,:),2);
            alla = [alla a];
            allb = [allb b];
            DATA.allskews(j,k,n) = skewness(b(:)-a(:));
        end
        DATA.skews(j,k) = skewness(allb(:)-alla(:));
        DATA.diffs{j,k} = cat(3,alla, allb);
    end
end
DATA.nsamples = nsmpls;
DATA.allx = allx;
figure(1);
imagesc(DATA.skews,'buttondownfcn',{@HitImage});
set(gca,'xtick',[1:length(DATA.nsamples)],'xticklabel',num2str(DATA.nsamples'));
xlabel('Nsamples group 1');
ylabel('Nsamples group 2');
set(gca,'ytick',[1:length(DATA.nsamples)],'yticklabel',num2str(DATA.nsamples'))
set(gcf,'UserData',DATA);
title(sprintf('Underlying distribution Skewness %.2f',skewness(allx)));
text(2,1,'Press here to redraw underlying distribtuion','buttondownfcn',{@HitImage},'color','w');
text(3,2,'Hit Squares to see resampled difference distribtuion','buttondownfcn',{@HitImage},'color','w');
text(4,3,'Shift Hit to see separate  distribtuions','buttondownfcn',{@HitImage},'color','w');
figure(2);
hist(allx,100);
title(sprintf('Underlying distribution Skewness %.2f',skewness(allx)));


function HitImage(a,b)
DATA = get(gcf,'UserData');
pos = get(gca,'currentpoint');
j = round(pos(1,2));
k = round(pos(1,1));
shiftalt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
figure(2);
hold off;
if k > j %hit in upper quadrant = show original distribution
    hist(DATA.allx(:),1000);
    title(sprintf('Underlying distribution Skewness %.2f',skewness(DATA.allx)));
else
    a = squeeze(DATA.diffs{j,k}(:,:,1));
    b = squeeze(DATA.diffs{j,k}(:,:,2));
    diffs = b-a;
    skew(1) = skewness(a(:));
    skew(2) = skewness(b(:));
    if shiftalt == 3 %shift hit
        xs = linspace(min([min(a(:)) min(b(:))]), max([max(a(:)) max(b(:))]));
        [y,x] = hist(a(:),xs);
        bar(x,y);
        [y,x] = hist(b(:),xs);
        y = y.* 1;
        hold on; plot(x,y,'r');
        xlabel('Resampled differences');
    elseif shiftalt == 2 %cntrl hit
            hist(DATA.allskews(j,k,:));
            title(sprintf('Samples %d,%d Skewness %.2f,%.2f. of difference %.2f',DATA.nsamples(j),DATA.nsamples(k),...
                skew(1),skew(2),skewness(diffs(:))));
            xlabel('skew in resample run for one sample set');
            hold on;
            x = mean(DATA.allskews(j,k,:));
            plot([x x],get(gca,'ylim'));
            x = skewness(diffs(:));
            plot([x x],get(gca,'ylim'),'r');
            text(x,mean(get(gca,'ylim')),'Skew of sum','color','r');
    else
        hist(diffs(:),100);
    end
skew(1) = skewness(a(:));
skew(2) = skewness(b(:));
title(sprintf('Samples %d,%d Skewness %.2f,%.2f. of difference %.2f',DATA.nsamples(j),DATA.nsamples(k),...
    skew(1),skew(2),skewness(diffs(:))));
end


function resampled = Bresample(values,nrep)
% Resample(values,nrep,nout)

n = length(values);
nout = n;

% Pick n random numbers between 1 and n: these will be the resampled values
newindices = unidrnd(n,nrep,nout);
resampled = values(newindices);
