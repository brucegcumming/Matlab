function [x, details] = TrackBlanks(Expt, varargin)
      
sid = [];
j = 1;
while j <= length(varargin)
    if strnmcpi(varargin{j},'tsamples',4)
        j = j+1;
        sid = varargin{j};
    end
    j = j+1;
end
blocks = Expt.Header.BlockStart;
tn = [Expt.Trials.Trial];
bExpt = Expt;
[nc,nr] = NsubPlots(length(blocks)-1);
nc = length(blocks)-1;
np = 1;
iprobe = [1:0.1:24];
res = PlotRevCorAny(Expt,'lfp');
[a, maxv] = max(res.sdfs.extras{1}.lfp,[],2);
[a, maxt] =  max(a);
if isempty(sid)
    sid = [maxt-5 maxt maxt+5];
end
for j = 1:length(Expt.Header.BlockStart)-1
    trials = find(tn >= blocks(j) & tn <= blocks(j+1));
    if length(trials) > 20
        subplot(1,nc,np);
        np = np+1;
        bExpt.Trials = Expt.Trials(trials);
        res = PlotRevCorAny(bExpt,'lfp');
        imagesc(res.sdfs.extras{1}.lfp');
        xtimes(j) =Expt.Trials(trials(1)).TrialStart./(10000 * 60);
        title(sprintf('ed %.2f, %.1fmin',Expt.Header.depths(j),xtimes(j)));
        a = res.sdfs.extras{1}.lfp(maxt,:);
        [a, maxs(j)] =  max(interp1([1:24],a,iprobe,'spline'));
        for k = 1:length(sid)
            zc = diff(sign(res.sdfs.extras{1}.lfp(sid(k),:)));
            id = find(zc(1:23) < 0);
            id = id(end);
            x(j,k) = interp1(res.sdfs.extras{1}.lfp(sid(k),id:id+1),[id id+1],0);
        end
    end
      end
  
  details.xtimes = xtimes;
 details.maxs = iprobe(maxs);