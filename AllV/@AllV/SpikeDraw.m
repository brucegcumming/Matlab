function SpikeDraw(a,b, mode)DATA = GetDataFromFig(a);DataClusters = AllV.mygetappdata(DATA,'Clusters');onoff = {'off' 'on'};redraw = 0;if mode == 2    DATA.spksperview = round(DATA.spksperview/2);elseif mode == 1    DATA.spksperview = DATA.spksperview *2;elseif mode == 3   DATA.plotspk.allprobes = 0;   AllV.SpoolSpikes(DATA);elseif mode == 4    DATA.plotdvdt = ~DATA.plotdvdt;    set(a,'Checked',onoff{DATA.plotdvdt+1});    AllV.PlotMeanSpike(DATA);elseif mode == 5    if DATA.plotcsd == 1        DATA.plotcsd = 0;    else        DATA.plotcsd = 1;    end    set(a,'Checked',onoff{DATA.plotcsd+1});    AllV.PlotMeanSpike(DATA);elseif strcmp(mode,'dvdy')    DATA.plotcsd =  ~(DATA.plotcsd == 2)*2;    set(a,'Checked',onoff{DATA.plotcsd/2+1});    redraw = 1;    AllV.PlotMeanSpike(DATA);elseif strcmp(mode,'subtrigger')    DATA.plotspk.subtrigger = ~DATA.plotspk.subtrigger;    set(a,'Checked',onoff{DATA.plotspk.subtrigger+1});    redraw = 1;    AllV.PlotMeanSpike(DATA);elseif mode == 6    DATA.plotspk.submean = ~DATA.plotspk.submean;    set(a,'Checked',onoff{DATA.plotspk.submean+1});    AllV.PlotMeanSpike(DATA);elseif mode == 7    DATA.plotspk.submax = ~DATA.plotspk.submax;    set(a,'Checked',onoff{DATA.plotspk.submax+1});    AllV.PlotMeanSpike(DATA);elseif mode == 8    DATA.plotspk.submin = ~DATA.plotspk.submin;    set(a,'Checked',onoff{DATA.plotspk.submin+1});    AllV.PlotMeanSpike(DATA);elseif mode == 9    DATA.plotspk.oneprobe = ~DATA.plotspk.oneprobe;    set(a,'Checked',onoff{DATA.plotspk.oneprobe+1});    AllV.PlotMeanSpike(DATA);elseif mode == 10    DATA.plotspk.bytrial = ~DATA.plotspk.bytrial;    set(a,'Checked',onoff{DATA.plotspk.bytrial+1});    AllV.SpoolSpikes(DATA);elseif mode == 11  %restrict time range    if ~isempty(DATA.restricttimerange)        t(1) = DATA.restricttimerange(1);    else        t(1) = DATA.t(1)-0.01;    end    t(2) = DATA.t(DATA.spklst(end))+0.01;    DATA = AllV.RestrictTimeRange(DATA,t);    DATA = AllV.SetPCs(DATA,1,0);elseif strcmp(mode,'spoolwithmean')    DATA.plotspk.showmean = ~DATA.plotspk.showmean;    set(a,'Checked',onoff{DATA.plotspk.showmean+1});    redraw = 1;elseif strcmp(mode,'allmeans')    AllV.PlotAllMeans(DATA);elseif strcmp(mode,'menu') % hit menu itself    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})elseif strcmp(mode,'excludetrials')    f = findobj('type','figure','Tag',DATA.tag.spikes);    it = findobj(f, 'Tag','ChooseTrial');    id = get(it,'value');    DATA = AllV.ExcludeTrials(DATA,id, 1);    DATA = AllV.SetPCs(DATA,1,0);elseif mode == 12 %restrict time range    if ~isempty(DATA.restricttimerange)        t(2)= DATA.restricttimerange(2);    else        t(2) = DATA.t(end)+0.01;    end    t(1)= DATA.t(DATA.spklst(1))-0.01;    DATA = AllV.RestrictTimeRange(DATA,t);    DATA = AllV.SetPCs(DATA,1,0);elseif mode == 13 %Spool All Probes    DATA.plotspk.allprobes = 1;    AllV.SpoolAllSpikes(DATA);elseif mode == 14 %Ues all spikes again    DATA = AllV.UseAllEvents(DATA);    set(DATA.toplevel,'UserData',DATA);elseif mode == 15 %Xcorr for selected probes    AllV.SetFigure(DATA.tag.xcorr, DATA,'front');    ps = find(DATA.selectprobe);    np = length(ps);    for j = 1:np        for k = 1:j            P = DataClusters{ps(j)};            Q = DataClusters{ps(k)};            xc = xcorrtimes(P.times,Q.times);            if k == j                xc(201) = 0;            end            subplot(np,np,k+(j-1)*length(ps));            bar(-200:200,xc);            axis('tight');            if k == j                title(sprintf('P%d',ps(j)));            end        end    endelseif strcmp(mode,'allquickspks') %QuickSpks All probes    DATA.plotspk.allprobes = 1;    AllV.SetFigure(DATA.tag.allspikes, DATA);    AllV.PlotQuickSpikes(DATA,100);elseif strcmp(mode,'allxy') %QuickSpks All probes    AllV.SetFigure(DATA.tag.allxy, DATA);    AllV.PlotAllProbes(DATA,'xy');elseif strcmp(mode,'xcorradj')  %Xcorr for all adjacent probes    xpts = linspace(DATA.spts(1),DATA.spts(end),401);    probes = SetProbesToUse(DATA, DATA.currentprobe, 2);    nc = 0;    for j = 1:length(probes)        nc = nc+1;        cells(nc).p = probes(j);        cells(nc).cl = 1;        for k = 1:length(DataClusters{cells(nc).p}.next)            nc = nc+1;            cells(nc).p = probes(j);            cells(nc).cl = k+1;        end    end    AllV.SetFigure(DATA.tag.xcorr,DATA);  AllV.PlotAllXCorr(DATA, DataClusters,cells);elseif strcmp(mode,'xcorrall') %Xcorr for all adjacent probes    nc = 1;    for k = 1:length(DataClusters)        if DATA.plot.dprimemin == 0 || (DataClusters{k}.fitdprime(1) < DATA.plot.dprimemin)        cells(nc).p = k;        cells(nc).cl = 1;        nc = nc+1;        end        for j = 1:length(DataClusters{k}.next)            if isfield(DataClusters{k}.next{j},'times') && ...                    (DATA.plot.dprimemin == 0 || DataClusters{k}.next{j}.fitdprime(1) < DATA.plot.dprimemin)                cells(nc).p = k;                cells(nc).cl = 1+j;                nc = nc+1;            end        end    end    AllV.SetFigure(DATA.tag.xcorr,DATA);    DATA.xcorrs = AllV.PlotAllXCorr(DATA, DataClusters, cells)endif ismember(mode,[1 2 4 5 6 7 8 9]) | redraw    DATA.spklst = DATA.spklst(1):DATA.spklst(1)+DATA.spksperview-1;    AllV.PlotSpikes(DATA,DATA.spklst);    set(DATA.toplevel,'UserData',DATA);end