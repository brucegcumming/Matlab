function DATA = CombineAutoCut(DATA,eid, expno, varargin)
 
     minimize = 0;
     plotmode = 1;
     show = 2;
     j = 1;
     mode = DATA.plot.autoclustermode;
     while j <= length(varargin)
         if strncmpi(varargin{j},'cluster',4)
             mode = 1;
         elseif strncmpi(varargin{j},'minimize',4)
             minimize = 1;
         elseif strncmpi(varargin{j},'noplot',4)
             plotmode = 0;
             show = 0;
         end
         j = j+1;
     end
     
     

     [DATA, ispk] = SetExptSpikes(DATA, eid, 'setrange');

     if length(ispk) < 4
         dprime = 0;
     end
     if length(ispk) > 3
         DATA = CalcClusterVars(DATA, DATA.Expts{eid}.gui.spkrange);

         if mode == 1
             if isfield(DATA,'AllClusters')
                 if iscell(DATA.AllClusters)
                     Cy = DATA.AllClusters{eid}(DATA.probe).cy;
                     Cx = DATA.AllClusters{eid}(DATA.probe).cx;
                 else
                 Cy = DATA.AllClusters(DATA.probe).cy;
                 Cx = DATA.AllClusters(DATA.probe).cx;
                 end
             else
                 Cy = DATA.Spikes.cy;
                 Cx = DATA.Spikes.cx;
             end
             [x,y,z] = CombineCalcDensity(DATA, ispk, 3);
             xv = x(1,:);
             yv = y(:,1);
             % find local maxima in smoothed 2-D density
             sd = 5;
             np = 10;

             nl = 1;
             stepscale = 1;
             while (np > 2 || np == 1) & stepscale > 0.1 & sd < 25 & sd > 0.25
                 if sd > 15
                     [a,b,K] = gauss2d(sd,[-60:60]);
                 else
                     [a,b,K] = gauss2d(sd,[-30:30]);
                 end
                 K = K./sum(K(:));
                 Z = conv2(z,K,'same');
                 dy = diff(Z,[],2);
                 dx = diff(Z);
                 yz = diff(sign(dy),[],2);
                 xz = diff(sign(dx));
                 lastnp = np;
                 pid = find(yz(2:199,:) < 0 & xz(1:198,2:199) < 0);
                 if length(ispk) > 5000
                     zthresh = max(Z(:))./40;
                 elseif length(pid) > 4
                     zthresh = max(Z(:))/10;
                 elseif length(pid) > 2
                     zthresh = max(Z(:))./20;
                 else
                     zthresh = 0;
                 end
                 pid = find(yz(2:199,:) < 0 & xz(1:198,2:199) < 0 & Z(2:199,2:199) > zthresh);

                 np = length(pid);
                 if np == 2
                     [py,  px] = ind2sub([198 198], pid);
                     d = abs((px(2)-px(1) + i * (py(2)-py(1))));
                     if d < 30
                         np = 1;
                     end
                 end
                 nps(nl) = np;
                 sds(nl) = sd;
                 zts(nl) = zthresh;
                 if nl ==1 && np > 10
                     stepscale = 4;
                 end
                 if np == 1 && sd < 50 && sd > 1
                     sd = sd - 0.5 * stepscale;
                 elseif np == 1 && sd  < 1
                     sd = sd * 0.8;
                 elseif np == 1
                     sd = sd * 0.8;
                 elseif nl > 1 && nps(nl-1) == 1
                     stepscale = stepscale  * 0.25;
                     sd = sd + stepscale;
                 else
                     sd = sd +stepscale;
                 end
                 nl = nl+1;
             end
             if np == 0
                 [a,pid] = max(Z(2:199,2:199));
                 np = 1;
             end
             [py,  px] = ind2sub([198 198], pid);
             px = px+1;
             py = py+1;
             [pd, peaki] = max(abs(px +i*py)); %find most distant local maximum
             peaka = [px(peaki) py(peaki)];

             % now find highest peak that is not this one
             [id, idx] = setdiff(pid,pid(peaki));
             zid = sub2ind([200 200],py(idx),px(idx));
             [a,b] = max(Z(zid));
             if np > 1
                 % This didn't help...
                 %             [a,c,K] = gauss2d(sd/2,[-30:30]);
                 %             K = K./sum(K(:));
                 %             Zf = conv2(z,K,'same');
                 peakb = [px(idx(b)) py(idx(b))];
                 xdi = linspace(x(1,peaka(1)),x(1,peakb(1)));
                 ydi = linspace(y(peaka(2),1),y(peakb(2),1));
                 dd = interp2(x,y,Z,xdi,ydi);
                 [dmin,b] = min(dd);
                 minpt = [xdi(b) ydi(b)];
                 zmax = Z(peaka(2),peaka(1));

                 d = sqrt(sum((peaka-peakb).^2));
             else
                 zmax = max(Z(:));
                 peakb = [];
                 d = 0;
             end
             if length(ispk) < 500
                 minzt = 0.1;
             elseif length(ispk) < 2000
                 minzt = 0.05;
             else
                 minzt = 0.02;
             end
             if np == 0
             elseif np == 1 || d < 30 || dmin/zmax > 0.8 || zmax/max(Z(:)) < minzt % too close treat as a single blob
                 if np > 1 %% peak is not for real
                     [a,id] = max(Z(:));
                     [peaka(2),  peaka(1)] = ind2sub([200 200], id);
                     if dmin/zmax > 0.8
                         dmin = zmax/5;
                     else
                         dmin = zmax/20;
                     end
                 else
                     dmin = zmax/50;
                 end
                 twopeaks = 0;
                 xu = max(Cx(ispk));
                 xl = peaka(1);
                 xc = (xv(xl)+xu)/2;
                 xr = (xu-xv(xl))/2;
                 id = find(Z(peaka(2):end,peaka(1)) < dmin);
                 if isempty(id)
                     id = find(Cx(ispk) > xv(xl));
                     yu = max(Cy(ispk(id)));
                 else
                     yu = yv(id(1)+peaka(2)-1);
                 end
                 yr = yu-yv(peaka(2));
                 yc = yv(peaka(2));
             else
                 twopeaks = 1;
                 %dmin is the low density value between teh cluster peak and the rest
                 %should find lower values in the outward directions, but may not - if the
                 %cluster is very good, dmin can be low.
                 id = find(Z(peaka(2),1:peaka(1)) > dmin);
                 if isempty(id)
                     id = find(Z(peaka(2),1:peaka(1)) < zmax/20);
                 end
                 if isempty(id)
                     xl = 1;
                 else
                     xl = id(1);
                 end
                 id = find(Z(peaka(2),peaka(1):end) < dmin);
                 if isempty(id)
                     xu = length(xv);
                 else
                     xu = id(1)+peaka(1)-1;
                 end
                 id = find(Z(1:peaka(2),peaka(1)) > dmin);
                 if isempty(id) % check for local minima
                     yl = 1;
                 else
                     yl = id(1);
                 end
                 id = find(Z(peaka(2):end,peaka(1)) < dmin);
                 if isempty(id)
                     yu = length(yv);
                 else
                     yu = id(1)+peaka(2)-1;
                 end
                 k = (xv(xu)-xv(xl))/(yv(yu)-yv(yl)); %ratio of radii
                 yr = sqrt((ydi(b)-ydi(1)).^2 + (xdi(b)-xdi(1)).^2/k^2);
                 xr = yr * k;

                 id = find(Z(peaka(2), peaka(1):end) < zmax/20);
                 if isempty(id)
                     xu = length(xv);
                 else
                     xu = peaka(1)+id(1);
                     if xu > length(xv)
                         xu = length(xv);
                     end
                 end
                 id = find(Z(peaka(2):end,peaka(1)) < zmax/20);
                 if isempty(id)
                     yu = max(Cy(ispk));
                 else
                     yu = yv(peaka(2)+id(1)-1);
                 end
                 k =  (xv(xu)-xdi(1))/xr;
                 xc = xdi(1) + xr * (k-1)/2;
                 xr = xr + xr *(k-1)/2;
                 k = (yu-ydi(1))/yr;
                 yc = ydi(1) + yr * (k-1)/2;
                 yr = yr + yr *(k-1)/2;
             end
             C.x = [xc xr 1];
             C.y = [yc yr 1];
             C.angle = 0;
             C.firstspk = NaN;
             C.lastspk = NaN;
             C.autocut = 2;
             C.params(1) = DATA.plot.clusterX;
             C.params(2) = DATA.plot.clusterY;
             C.Arange = DATA.clusterArange;
             C.Brange = DATA.clusterBrange;
             C.Erange = DATA.clusterErange;
             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA.Expts{eid}.gui.classified = 3;
             if plotmode
                 GetFigure('AutoCut');
                 subplot(2,1,1);
                 hold off;
                 pcolor(x,y,Z); shading('interp');
                 hold on;
                 plot(xv(peaka(1)),yv(peaka(2)),'ok');
                 if ~isempty(peakb)
                     plot(xv(peakb(1)),yv(peakb(2)),'ow');
                 end
                 plot(xc,yc,'+k');
                 if np > 1
                     plot(minpt(1),minpt(2),'+w');
                 end
                 x = linspace(0,xr);
                 y =  sqrt(yr.^2 - (x.*yr/xr).^2);
                 x = [x fliplr(x) -x fliplr(-x)];
                 y = [y fliplr(-y) -y fliplr(y)];
                 plot(x+xc,y+yc,'r');
                 title(sprintf('%s%d:P%d %d peaks',DATA.Expts{eid}.Header.expname,eid,DATA.probe,np));
                 subplot(2,1,2);
                 hold off;
                 plot(Z(peaka(1),:));
                 hold on;
                 plot(Z(:,peaka(2)),'r');
                 if isfigure(DATA.figs.allprobes)
                     figure(DATA.figs.allprobes);
                     [nr,nc] = Nsubplots(length(DATA.probelist));
                     subplot(nr,nc,DATA.probe);
                 else
                     GetFigure(DATA.tag.clusterxy);
                 end
                 [DATA, ispk, dprime, details] = SetExptSpikes(DATA, eid, show,'useexpall');
             else
                 [DATA, ispk, dprime, details] = SetExptSpikes(DATA, eid, 0,'useexpall');

             end
             DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
             nspk = sum(details.nc);

             % IF there are two peaks, optimize cluster size
             if twopeaks & nspk > 0
                 oldC = C;
                 if minimize
                     DATA.ispk = ispk;
                     C = OptimizeDprime(DATA);
                     DATA.cluster{1,DATA.probe} = C;
                     [DATA, newd, d] = SetSpkCodes(DATA,ispk,DATA.probe,0);
                 else
                     newd = dprime * 1.001;
                     dps(1) = dprime;
                     nl = 1;
                     frac = 0;
                     while newd > dprime
                         dprime = newd;
                         xr = xr.*1.1;
                         %                 yr = yr.*1.1;
                         C.x = [xc xr 1];
                         C.y = [yc yr 1];
                         DATA.cluster{1,DATA.probe} = C;
                         [DATA, newd, d] = SetSpkCodes(DATA,ispk,DATA.probe, 0);
                         frac = d.nc(1)./length(ispk);
                         nl = nl+1;
                         dps(nl) = newd;
                     end
                     C.x = [xc xr./1.1 1];
                     dprime = dps(nl-1);
                     newd = dprime*1.001;
                     while newd > dprime
                         dprime = newd;
                         %                 xr = xr.*1.1;
                         yr = yr.*1.1;
                         %                C.x = [xc xr 1];
                         C.y = [yc yr 1];
                         DATA.cluster{1,DATA.probe} = C;
                         [DATA, newd, d] = SetSpkCodes(DATA,ispk,DATA.probe,0);
                         nl = nl+1;
                         dps(nl) = newd;
                     end
                     C.y = [yc yr./1.1 1];
                 end
                 frac = d.nc(1)./length(ispk);
                 if frac > 0.99  %didn't work
                     C = oldC;
                 end

                 DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
                 DATA.cluster{1,DATA.probe} = C;
                 [DATA, dprime, details] = SetSpkCodes(DATA,ispk,DATA.probe,show);
                 DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
             end
             if np == 1 | dprime < 1
                 C.dprime = dprime;
                 DATA = CheckAutoRate(DATA, eid, C, 10);
             end
             if show
%                 DrawClusters(DATA,DATA.Expts{eid}.Cluster,0);
             end
             title(sprintf('%s%d %d/%d spikes: dprime %.1f',DATA.Expts{eid}.Header.expname,eid,details.nc(1),length(ispk),dprime));

             drawnow;
             return;





         elseif mode == 2
             [a,b] = max(z);
             py = smooth(y(b,1),5,'gauss');
             dens = smooth(a,5,'gauss');
             [a,maxi] = max(dens);
             maxpt = [b(maxi) maxi];
             for j = 1:length(dens)
                 dip(j) = dens(j)./max(dens(1:j));
             end
             id = find(diff(sign(diff(dip))) >0);
             GetFigure('AutoCut');
             subplot(2,1,1);
             hold off;
             plot(x(1,:),dip);
             hold on;
             for j = 1:length(id)
                 plot([x(1,id(j)+1) x(1,id(j)+1)],[0 1],':');
             end
             id = find(dip(2:end) > 0.5 & diff(dip) > 0.005);
             [a,b] = min(dip(id));
             if isempty(id)
                 id = find(dip < max(dip)/2);
                 ll = id(1);
                 ul = length(dens);
             else
                 ll = id(b);
                 %     ll = id(end);
                 ul= id(end)+10;
                 id = find(diff(dip(1:ll)) <= 0);


                 ll = id(end); %first positive ddip
                 peakdip = max(dip(ll:end))
                 id = find(diff(dip(ul:end)) > 0 & dip(ul+1:end) > peakdip);
                 if id
                     ul = ul+id(1);
                 else
                     ul = length(dens);
                 end
             end

             if x(1,ul) > max(Cx(ispk)) | ul == size(x,1)
                 xu = max(Cx(ispk)).*1.05;
             else
                 xu = x(1,ul);
             end
             plot([x(1,ll) x(1,ll)],[0 1],'r:');
             plot([xu xu],[0 1],'r:');
             subplot(2,1,2);
             hold off;
             pcolor(x,y,z);
             hold on;
             shading('flat');
             plot(x(1,:),py);

             xc = mean([xu x(1,ll)]);
             xr = (xu-x(1,ll))/2;
             yd = sum(z(:,ll:ul),2);
             yc = sum(y(:,1).*yd)./sum(yd);
             sigma = sqrt(sum(yd.*(y(:,1)-yc).^2)./sum(yd));
             yr = 6.* sigma;
             id = find(Cx > x(1,ll) & Cx <= xu);
             sigma = std(Cy(id));

             C.x = [xc xr 1];
             C.y = [yc yr 1];
             C.angle = 0;
             C.firstspk = NaN;
             C.lastspk = NaN;
             C.autocut = 2;
             C.params(1) = DATA.plot.clusterX;
             C.params(2) = DATA.plot.clusterY;
             C.Arange = DATA.clusterArange;
             C.Brange = DATA.clusterBrange;
             C.Erange = DATA.clusterErange;

             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA.Expts{eid}.gui.classified = 3;
             DATA.Expts{eid}.gui.clustertype = 3;

             GetFigure(DATA.tag.clusterxy);
             [DATA, ispk, dprime]  = SetExptSpikes(DATA, eid, show,'useexpall');
             if np == 1 | dprime < 2
                 DATA = CheckAutoRate(DATA, eid, C, 10);
             end
             DrawClusters(DATA,DATA.Expts{eid}.Cluster,0);
             return;
         end
         xid = find(DATA.Spikes.cx(ispk) > median(DATA.Spikes.cx(ispk)));
         yid = find(DATA.Spikes.cy(ispk) < median(DATA.Spikes.cy(ispk)));
         xc = prctile(DATA.Spikes.cx(ispk(yid)),90);
         xmax = prctile(DATA.Spikes.cx(ispk(yid)),99);
         yc = prctile(DATA.Spikes.cy(ispk(xid)),90);
         ymax = prctile(DATA.Spikes.cy(ispk(xid)),99);
         xr = (xmax-xc) * 2;
         if xr == 0
             xr = xc/3;
         end
         yr = (ymax-yc) * 1.5;
         if yr == 0
             yr = yc/3;
         end
         xc = xc.*1.4;
         yc = yc.*0.7;
         C.x = [xc xr 1];
         C.y = [yc yr 1];
         C.angle = 0;
         C.firstspk = NaN;
         C.lastspk = NaN;
         C.autocut = 1;
         DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
         if show & isfigure(DATA.xyfig);
             GetFigure(DATA.xyfig);
             hold off;
         end
         if ~isfield(DATA,'ptsize')
             DATA.ptsize = 4;
         end
         DATA = SetExptSpikes(DATA, eid, show,'useexpall');
         dprime = DATA.cluster{1,DATA.probe}.dprime;
         [rate, frac] = SumExptSpikes(DATA, eid, 1);
         if isinf(xr) | isinf(yr)
             fprintf('xc,yc Inf problem %f %f',xc,yc);
             return;
         end
         while rate < 10 & frac < 0.9
             xr = xr .*1.2;
             yr = yr.*1.2;
             if ~(xr < 1000) | ~(yr < 1000) %catches inf also..
                 fprintf('Big raduis',xc,yc);
             end
             C.x = [xc xr 1];
             C.y = [yc yr 1];
             fprintf('Only %.2f Hz (%.2f) Trying %.2f,%.2f\n',rate,frac,xr,yr);
             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA = SetExptSpikes(DATA, eid, show,'useexpall');
             dprime = DATA.cluster{1,DATA.probe}.dprime;
             [rate, frac] = SumExptSpikes(DATA, eid, 1);
         end
         while rate > 100 & frac > 0.1
             xc = xc .* 1.2;
             yc = yc .* 1.2;
             C.x = [xc xr 1];
             C.y = [yc yr 1];
             fprintf('Rate %.2f Hz (%.2f) Trying %.2f,%.2f\n',rate,frac,xc,yc);
             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA = SetExptSpikes(DATA, eid, show,'useexpall');
             dprime = DATA.cluster{1,DATA.probe}.dprime;
             [rate, frac] = SumExptSpikes(DATA, eid, 1);
         end
         DATA.cluster{1,DATA.probe}.autocut = 1;
         DATA.cluster{1,DATA.probe}.dprime = dprime;
         DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
         fprintf('Probe %d Expt %d Cluster at %.3f+- %.3f, %.3f +- %.3f %.2fHz\n',DATA.probe,expno, C.x(1), C.x(2), C.y(1), C.y(2),rate);
         if show
             title(sprintf('probe %d, Expt %d %.f Hz Dp %.2f',DATA.probe,expno,rate,dprime));
         end
     else
         fprintf('Probe %d Expt %d No Spikes\n',DATA.probe,expno);
         C.firstspk = NaN;
         C.lastspk = NaN;
         C.autocut = 1;
         DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
     end

     function DATA = CheckAutoRate(DATA, eid, C, minrate)
     xr = C.x(2);
     yr = C.y(2);
     xc = C.x(1);
     yc = C.y(1);
     show = 2;
     [rate, frac, nspk] = SumExptSpikes(DATA, eid, 1);
     ispk = DATA.Expts{eid}.gui.spks;
     xm = max(DATA.Spikes.cx(ispk));
     while rate < minrate & frac < 0.9
         xr = xr .*1.2;
         yr = yr.*1.2;
         if ~(xr < 1000) | ~(yr < 1000) %catches inf also..
             fprintf('Big raduis',xc,yc);
             rate = NaN;
         end
         C.x = [xc xr 1];
         C.y = [yc yr 1];
         fprintf('Only %.2f Hz (%.2f) Trying %.2f,%.2f\n',rate,frac,xr,yr);
         DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
         [DATA, ispk, dprime]  = SetExptSpikes(DATA, eid, show,'useexpall');
         DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
         DATA.Expts{eid}.Cluster{1,DATA.probe}.nspk = nspk;
         if ~isnan(rate)
         [rate, frac] = SumExptSpikes(DATA, eid, 1);
         end
     end
     
     
    function [meanrate, proportion, totalspks] = SumExptSpikes(DATA, expid, cluster)
    cspks = DATA.Expts{expid}.gui.spks;
    ctype = 2;
    j = 1;
    for trial = [DATA.Expts{expid}.Trials]
        if isfield(DATA,'AllClusters')
            if iscell(DATA.AllClusters)
                ispk = find(DATA.AllClusters{expid}(DATA.probe).times(cspks) > trial.Start(1) & DATA.AllClusters{expid}(DATA.probe).times(cspks) < trial.End(end)+500);
                ispks(j) = sum(DATA.AllClusters{expid}(DATA.probe).codes(cspks(ispk),1) == cluster);
            else
                ispk = find(DATA.AllClusters(DATA.probe).times(cspks) > trial.Start(1) & DATA.AllClusters(DATA.probe).times(cspks) < trial.End(end)+500);
                ispks(j) = sum(DATA.AllClusters(DATA.probe).codes(cspks(ispk)) == cluster);
            end
        elseif isfield(DATA,'AllSpikes')
        ispk = find(DATA.AllSpikes{DATA.probe}.times(cspks) > trial.Start(1) & DATA.AllSpikes{DATA.probe}.times(cspks) < trial.End(end)+500);
        ispks(j) = sum(DATA.AllSpikes{DATA.probe}.codes(cspks(ispk),ctype) == cluster);
        else
        ispk = find(DATA.AllData.Spikes.times(cspks) > trial.Start(1) & DATA.AllData.Spikes.times(cspks) < trial.End(end)+500);
    ispks(j) = sum(DATA.AllData.Spikes.codes(cspks(ispk),ctype) == cluster);
        end
    aspks(j) = length(ispk);
    dur(j) = trial.End(end)-trial.Start(1);
    j = j+1;
    end
    meanrate = 10000.* mean(ispks)./mean(dur);
    proportion = sum(ispks)./sum(aspks);
    totalspks = sum(ispks);
  