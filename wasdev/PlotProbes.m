function PlotProbes(rcs, varargin)

plottype = 'LFPpwrdiff'
probesep = 100;
cptype = 1;
showtypes = [0 1 0 0];
plotbands = [3:4];
freqbands = [2 10 30 60 90];
figb = 'LFP power by choice';

j = 1;
while  j <= length(varargin)
    if strncmpi(varargin{j},'bands',4)
        j = j+1;
        plotbands = varargin{j};
    elseif strncmpi(varargin{j},'figb',4)
        j = j+1;
        figb = varargin{j};
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'types',4)
        j = j+1;
        showtypes = varargin{j};
    end

j = j+1;
end
plots = find(showtypes);

    for j = 1:length(freqbands)-1
    flabels{j} = sprintf('%d-%d',freqbands(j),freqbands(j+1));
end

         GetFigure('Psych');
            cellcp = []; cellpcp = []; celllcp = [];
            cellpos = [];
            cellcv = [];
            probes = [];
            cps = [];
            alln = zeros(100,1);
            allbands = zeros(4,100,4);
            for j = 1:length(rcs)
                pk(j,:) = rcs{j}.pk.kernel;
                pkx(j,:) = rcs{j}.pk.xv;
                upstim(j) = mean([rcs{j}.cp.dnstim]);
                probes(j,:) = (rcs{j}.probes - rcs{j}.blankmax);
                cps(j,:) = [rcs{j}.cp.cp];
                pcps(j,:) = [rcs{j}.cp.pcp];
                cvs(j,:) = rcs{j}.cv;
                allcp(j,:) = [rcs{j}.cp.cp];
                if isfield(rcs{j},'cellcp')
                cellcp = [cellcp rcs{j}.cellcp];
                celllcp = [celllcp rcs{j}.celllcp];
                cellpcp = [cellpcp rcs{j}.cellpcp];
                cellpos = [cellpos (rcs{j}.cellprobe - rcs{j}.blankmax).* probesep];
                cellcv = [cellcv rcs{j}.cellcv];
                rid(length(cellcp)) = 0;
                rid(find(rid == 0)) = j;
                end
                if strcmp(plottype,'LFPpwrdiff')
                    iprobes(j,:) = 20 +rcs{j}.probes - rcs{j}.blankmax;
                    allbands(:, iprobes(j,:),1,j) = rcs{j}.pwrbands;
                    s =  repmat(sign(sum(rcs{j}.pwrdiff)),4,1);
                    allbands(:, iprobes(j,:),2,j) = rcs{j}.pwrbands .*s;
                    s = abs(sin(rcs{j}.x(rcs{j}.prefs,1) .*pi/180));
                    s = s > sin(pi/4);
                    s = repmat(sign(s' -0.5),4,1);
                    allbands(:, iprobes(j,:),3,j) = rcs{j}.pwrbands .*s;
                    s = abs(sin(rcs{j}.prefdir .*pi/180));
                    s = s > sin(pi/4);
                    s = repmat(sign(s -0.5),4,1);
                    allbands(:, iprobes(j,:),4,j) = rcs{j}.pwrbands .*s;
                    alln(iprobes(j,:),:) = alln(iprobes(j,:),:) + 1;
                elseif strcmp(plottype,'sum')
                    iprobes(j,:) = 20 +rcs{j}.probes - rcs{j}.blankmax;
                    pwr = squeeze(sum(rcs{j}.lfpsigpwr));
                    pwr = squeeze(sum(rcs{j}.lfpwr));
                    for k = 1:length(freqbands)-1
                        id = find(rcs{j}.lfpfrq >= freqbands(k) & rcs{j}.lfpfrq < freqbands(k+1));
                            pwrbands(k,:) = sum(pwr(id,:))./length(id);
                            
                    end

                    allbands(:, iprobes(j,:),1,j) = pwrbands;
                    alln(iprobes(j,:),:) = alln(iprobes(j,:),:) + 1 ;
                end
            end
            subplot(3,1,1);
            hold off;
            [a,b] = sort(upstim);
            imagesc([min(pkx(:)) max(pkx(:))], 1:length(b), pk(b,:));
            hold on;
            plot(upstim(b),1:length(b),'w:');
            colorbar;
            subplot(3,1,3);
            hold off;
            muid = find(cvs > 0.2);
            id = find(cellcv > 0.7);
            if cptype == 1 %by pref
            plot(probes(muid), pcps(muid),'o');
            hold on;
            plot(cellpos(id),cellpcp(id),'ko','markerfacecolor','k');
            title(sprintf('Mean CP %.3f, Cells %.3f',mean(pcps(muid)),mean(cellpcp(id))));
            elseif cptype ==2
            plot(probes(muid), cps(muid),'o');
            hold on;
            plot(cellpos(id),cellcp(id),'ko','markerfacecolor','k');
            title(sprintf('Mean CP %.3f, Cells %.3f',mean(cps(muid)),mean(celllcp(id))));
            else
            plot(probes(muid), cps(muid),'o');
            hold on;
            plot(cellpos(id),cellcp(id),'ko','markerfacecolor','k');
            title(sprintf('Mean CP %.3f, Cells %.3f',mean(cps(muid)),mean(cellcp(id))));
            end
            
            plot([min(probes(muid)) max(probes(muid))],[0.5 0.5],'k:');
            if strcmp(plottype,'LFPpwrdiff')
                GetFigure(figb);
                id = find(alln > 5);
                offset = id(1) - 20;
                nm = repmat(alln(id)',[length(plotbands) 1 4]);
                allsd = squeeze(std(allbands(plotbands,id,:,:),[],4))./sqrt(nm);
                mband = squeeze(sum(allbands(plotbands,id,:,:),4))./nm;
                depths = ([1:size(mband,2)] +offset) .* probesep;
                depths = repmat(depths',1,size(mband,1));
                for j = 1:length(plots)
                subplot(length(plots),1,j);
                errorbar(depths,squeeze(mband(:,:,plots(j)))',squeeze(allsd(:,:,plots(j)))');
                if j == 1
                    legend(flabels{plotbands});
                end
            end
            elseif strcmp(plottype,'sum')
                GetFigure(figb);
                id = find(alln > 5);
                offset = id(1) - 20;
                nm = repmat(alln(id)',[length(plotbands) 1]);
                allsd = squeeze(std(allbands(plotbands,id,1,:),[],4))./sqrt(nm);
                mband = squeeze(sum(allbands(plotbands,id,1,:),4))./nm;
                depths = ([1:size(mband,2)] +offset) .* probesep;
                depths = repmat(depths',1,size(mband,1));
                errorbar(depths,mband',allsd');
        end
   