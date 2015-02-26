function result= PlotRCLFP(rc, Expt, varargin)
% Plot rc struct returned by PlotRevCorAny, for LFP
% ...,'chan', chan)   sets the LFP probe #, otherwise plots 1
%

%TO DO:
%subtract mean overall LFP to remove (a) framerate, (b) onset transient
%get peak time for spike var for the bystim analysis

trigbystim = 0; % set this to 1 to test Spike triggered lfp in viciniyt of stim
colors = mycolors;
lfptag = 'LFPRC';
j = 1;
submean = 0;
nopsych = 0;
lfpch = 1;
mkspktrig = 1;
yval = 1;
legendpos = 'Best';

while j <= length(varargin)
    if strncmpi(varargin{j},'bystim',4)
        trigbystim = 1;
    elseif strncmpi(varargin{j},'chan',3)
        j = j+1;
        lfpch = varargin{j};
    elseif strncmpi(varargin{j},'legendpop',3)
        j = j+1;
        legendpos = varargin{j};
    elseif strncmpi(varargin{j},'submean',4)
        submean = 1;
    elseif strncmpi(varargin{j},'yval',4)
        j = j+1;
        yval = varargin{j};
    end
    j = j+1;
end
GetFigure(lfptag);

if nargin < 2 | ~isfield(Expt.Trials,'Spikes')
    mkspktrig = 0;
end
if ~isfield(rc,'x') & isfield(rc.sdfs,'x')
    rc.x = rc.sdfs.x;
    rc.y = rc.sdfs.y;
end
subplot(2,1,1);
hold off;
for j = 1:size(rc.sdfs.lfp,1)
    for k = 1:size(rc.sdfs.lfp,2)
        sdfs(:,j,k) = rc.sdfs.lfp{j,k}(:,lfpch);
    end
end

meansdf = squeeze(mean(mean(sdfs,3),2)); %overall mean
if submean
    refsdf = meansdf;
else
    refsdf = zeros(size(meansdf));
end
result.meansdf = meansdf;
for j = 1:size(rc.sdfs.lfp,1)
    h(j) = plot(rc.sdfs.lfptimes./10,rc.sdfs.lfp{j,1}(:,lfpch)-refsdf,'color',colors{j});
    hold on;
    for k = 2:size(rc.sdfs.lfp,2)
    plot(rc.sdfs.lfptimes./10,rc.sdfs.lfp{j,k}(:,lfpch)-refsdf,'color',colors{j},'linestyle',':');
    end
    labels{j} = sprintf('%.1f',rc.x(j,yval));
    result.lfpall(j,:) = rc.sdfs.lfp{j,yval}(:,lfpch)-refsdf;
end
nl = j;
for j = 1:size(rc.sdfs.extras,1)
    h(j+nl) = plot(rc.sdfs.lfptimes./10,rc.sdfs.extras{j,1}.lfp(:,lfpch)-refsdf,'color',colors{j},'linestyle','--');
    hold on;
    for k = 2:size(rc.sdfs.extras,2)
    plot(rc.sdfs.lfptimes./10,rc.sdfs.extras{j,k}.lfp(:,lfpch)-refsdf,'color',colors{j},'linestyle','-.');
    end
    labels{j+nl} = sprintf('%.1f',rc.sdfs.extras{j,1}.label);
    result.extralfp(j,:) = rc.sdfs.extras{j,yval}.lfp(:,lfpch)-refsdf;
end
if ~submean
    plot(rc.sdfs.lfptimes./10,squeeze(mean(mean(sdfs,3),2)),'k'); %overall mean
end
[a,b] = max(var(sum(sdfs,3)'));
tmax = rc.sdfs.lfptimes(b);
idx = b-8:b+8;
legend(h,labels,'Location',legendpos);
subplot(2,1,2);
hold off;
%plot mean of lfp signal over range of times b
plot(rc.x(:,1),mean(sdfs(idx,:,1)));
hold on;
if size(sdfs,3) > 1
plot(rc.x(:,1),mean(sdfs(idx,:,2)),'r'); %Down Choices
end
spks = squeeze(mean(rc.y,2));
[a, tmax] = max(var(spks));
ymax = get(gca,'Ylim');
spscale = diff(ymax)./max(spks(:,tmax));
plot(rc.x(:,1),ymax(1)+rc.y(:,1,tmax).*spscale,':');
if size(sdfs,3) > 1
plot(rc.x(:,1),ymax(1)+rc.y(:,2,tmax).*spscale,'r:');
end



if isfield(Expt.Trials,'Dc')
for j = 1:length(Expt.Trials)
    Expt.Trials(j).Dc = Expt.Trials(j).Dc(end);
    if isfield(Expt.Trials,'Pd')
        Expt.Trials(j).sigPd = Expt.Trials(j).Pd(end);
    end
    lfprms(j,:) = sqrt(mean(Expt.Trials(j).LFP.^2));
end
maxlpwr = prctile(lfprms,90) * 3;
%
% for disparity, near (negative) is RespDir +1, = down choices
uid = find([Expt.Trials.RespDir] < 0 & [Expt.Trials.Dc] == 0 & lfprms(lfpch) < maxlpwr(lfpch));

ulfpvals = [];
if isfield(Expt.Trials,'Spikes')
for j = uid
    starts = Expt.Trials(j).lfpo+Expt.Trials(j).Spikes./(Expt.Header.LFPsamplerate * 10000);
    starts = round(starts);
    starts = starts(find(starts > 0 & starts < length(Expt.Trials(j).LFP)));
    ulfpvals = [ulfpvals; Expt.Trials(j).LFP(starts)];
end
did = find([Expt.Trials.RespDir] > 0 & [Expt.Trials.Dc] == 0 & lfprms < maxlpwr);
else
    nopsych = 1;
    plothist = 6;
end
end
dlfpvals = [];
delay = 950;
hw = rc.frameperiod/2;
type = 'or';

if mkspktrig
    if trigbystim
    %[result.aavg, result.ares] = SpTrigLFP(Expt.Trials,20000,1/(Expt.Header.LFPsamplerate * 10000),100);
    nstim = 1;
        for j = 1:length(Expt.Trials)
            spks = [];
            sid = find(ismember(Expt.Trials(j).(type),[0 30]));
            for k = sid';
                s = Expt.Trials(j).Start(k) - Expt.Trials(j).Start(1) + delay;
                spks = [spks find(Expt.Trials(j).Spikes' > s & ...
                    Expt.Trials(j).Spikes' < s + 106)];
            end
            nstim = nstim + length(sid);
            Expt.Trials(j).spkuse = spks;
            Expt.Trials(j).stdelays = Expt.Trials(j).Start(sid) - Expt.Trials(j).Start(1) + delay + hw;
        end
    [result.davg, dres] = SpTrigLFP(Expt.Trials,20000,1/(Expt.Header.LFPsamplerate * 10000),100,...
       'frameper',rc.frameperiod);
    dres.nstim = nstim;
    nstim = 1;
        for j = 1:length(Expt.Trials)
            spks = [];
            sid = find(ismember(Expt.Trials(j).(type),[90 120]));
            for k = sid';
                s = Expt.Trials(j).Start(k) - Expt.Trials(j).Start(1) + delay;
                spks = [spks find(Expt.Trials(j).Spikes' > s & ...
                    Expt.Trials(j).Spikes' < s + 106)];
            end
            Expt.Trials(j).spkuse = spks;
            Expt.Trials(j).stdelays = Expt.Trials(j).Start(sid) - Expt.Trials(j).Start(1) + delay + hw;
            nstim = nstim + length(sid);
        end
    [result.uavg, ures] = SpTrigLFP(Expt.Trials,20000,1/(Expt.Header.LFPsamplerate * 10000),...
        100,'histvals',dres.vx,'frameper',rc.frameperiod);
        ures.nstim = nstim;
    elseif nopsych
        [result.davg, dres] = SpTrigLFP(Expt.Trials,20000,1/(Expt.Header.LFPsamplerate * 10000),100);
        did = 1:length(Expt.Trials);
        plothist = 6;
        ures = [];
    else
    [result.davg, dres] = SpTrigLFP(Expt.Trials(did),20000,1/(Expt.Header.LFPsamplerate * 10000),100);
    [result.uavg, ures] = SpTrigLFP(Expt.Trials(uid),20000,1/(Expt.Header.LFPsamplerate * 10000),...
        100,'histvals',dres.vx);
    plothist = 5;
    end
result.dtrig = dres;
result.utrig =ures;

for j = did
    starts = Expt.Trials(j).lfpo+Expt.Trials(j).Spikes./(Expt.Header.LFPsamplerate * 10000);
        starts = round(starts);
    starts = starts(find(starts > 0 & starts < length(Expt.Trials(j).LFP)));
    dlfpvals = [dlfpvals; Expt.Trials(j).LFP(starts)];
end
else
    plothist = 6;
end
GetFigure('LFPhist');
subplot(2,1,1);
hold off;



if plothist == 1
[y,x] = smhist(ulfpvals);
plot(x,y./sum(y));
[ay,ax] = smhist(dlfpvals,'nsmp',length(y));
hold on;
plot(ax,ay./sum(ay),'r');
elseif plothist == 3
    [a,b] = sort(ures.framev);
    v = smooth(a,100);
    p = smooth(ures.counts(b),100);
    hold off;
    plot(v,p);
    [a,b] = sort(dres.framev);
    v = smooth(a,100);
    p = smooth(dres.counts(b),100);
    hold on;
    plot(v,p,'r');
elseif plothist == 4  %% plot F dists for V, for each frame spikecount
    % this is best way to see relationship. 
    [y,x] = smhist(dres.framev(find(dres.counts ==0)));
    hold off;
    plot(x,y);
    hold on;
    x = mean(dres.framev(find(dres.counts ==0)));
    plot([x x],[0 max(y)*1.1]);
    [y,x] = smhist(dres.framev(find(dres.counts ==1)));
    plot(x,y,'r');
    x = mean(dres.framev(find(dres.counts ==1)));
    plot([x x],[0 max(y)*1.1],'r');
    [y,x] = smhist(dres.framev(find(dres.counts ==2)));
    plot(x,y,'g');
    x = mean(dres.framev(find(dres.counts ==2)));
    plot([x x],[0 max(y)*1.1],'g');

    [y,x] = smhist(ures.framev(find(ures.counts ==2)));
    plot(x,y,'g--');
    x = mean(ures.framev(find(ures.counts ==2)));
    plot([x x],[0 max(y)*1.1],'g--');
    
    [y,x] = smhist(ures.framev(find(ures.counts ==1)));
    plot(x,y,'r--');
    x = mean(ures.framev(find(ures.counts ==1)));
    plot([x x],[0 max(y)*1.1],'r--');
    [y,x] = smhist(ures.framev(find(ures.counts ==0)));
    plot(x,y,'--');
    x = mean(ures.framev(find(ures.counts ==0)));
    plot([x x],[0 max(y)*1.1],'--');

    legend('0','1','2');
elseif plothist == 5  %% Plot firing vs V, for stim means + within
    hold off;
    shift = 10;
    scale = rc.frameperiod./10000;  % convert rate to spikes/frame
    for j = 1:length(rc.sdfs.lfp)
    y = interp1(rc.sdfs.lfptimes,rc.sdfs.lfp{j},rc.times);
    plot(rc.sdfs.s{j}(10:end-shift) .* scale, y(10+shift:end),'color',colors{j});
    hold on;
    end
    plot([0:length(dres.mv)-1],dres.mv,'o-');
    plot([0:length(ures.mv)-1],ures.mv,'ro-');
elseif plothist == 6 % no psych
    ures = [];
else
    plot(ures.vx,ures.vy./sum(ures.vy),':');
    hold on;
    plot(ures.vx,ures.trigvy./sum(ures.trigvy));
    plot(dres.vx,dres.vy./sum(dres.vy),'r:');
    plot(dres.vx,dres.trigvy./sum(dres.trigvy),'r');
end
title('LFP histogram');
legend('All','triggered');
subplot(2,1,2);
hold off;
if isfield(ures,'times')
plot(ures.times./10,result.uavg);
hold on;
end
if mkspktrig
plot(dres.times./10,result.davg,'r');
legend('Up','Down');
title('Spike Triggered');
if isfield(result,'uavg')
[a,b,ai,bi] = fit_bothsubj2error(result.davg,result.uavg);
end
end