function CompareSpikeWaves(P,Q,aid,bid, varargin)
%CompareSpikeWaves(P,Q,aid,bid) Compare a set of spikes 
%from two spks structs P,!, and two vectors of indices.
%Currenlty length(aid) mu s= length(bid)
%used e.g. by for PC.syncspikes
%...,'clusterid',clnum) give cluter #s of spikes
voff = P.maxv./2;
xpos = 0;
x = [1:size(P.values,2)]./40;
labels{1} = sprintf('P%d',P.probe);
labels{2} = sprintf('P%d',Q.probe);
for j = 1:length(varargin)
    if strncmpi(varargin{j},'clusterid',6)
        j = j+1;
        clnum = varargin{j};
        labels{1} = sprintf('P%d/%d',P.probe,clnum(1));
        labels{2} = sprintf('P%d/%d',Q.probe,clnum(2));
    elseif strncmpi(varargin{j},'xpos',4)
        j = j+1;
        xpos = varargin{j};
        x = x+xpos;
    end
    j = j+1;
end
if length(aid) ~= length(bid) && 0
    return;
end
xp = 0;
yp = 0;
Q = Spk2Double(Q);
P = Spk2Double(P);
if isfield(Q,'xchans') %find overlapping xchans
         [a,b] = ismember(Q.xchans,P.xchans);
         xp = find(ismember(Q.xchans,P.probe));
         if isempty(xp)
             xp = find(b);
             [a,b] = ismember(P.xchans,Q.xchans);
             yp = find(b);
             voff = voff;
         end
         if isempty(xp)
             xp = 0;
         elseif length(xp) > 1
             xp = xp(1);
             yp = yp(1)
         end
end

if length(aid) == length(bid)
    toff = (P.times(aid)-Q.times(bid))./10;
else
    toff = 0;
end
plot(x,P.values(aid,:),'k');
hold on;
X = bsxfun(@minus,x,toff);
plot(X',double(Q.values(bid,:))'+voff,'b');
if isfield(Q,'xchans')
    if xp && yp %a shared extra channel
        vo = max([max(max(Q.xvalues(xp,:,bid))) max(max(P.xvalues(yp,:,aid))) -min(min(P.values(aid,:)))]);
        plot(X',squeeze(double(Q.xvalues(xp,:,bid)))-vo.*0.98,'b');
        plot(x,squeeze(double(P.xvalues(yp,:,aid)))-vo,'k');
        text(x(end),-vo,sprintf('P%d',Q.xchans(xp)),'horizontalalignment','left','color','b','fontsize',20);
    elseif xp %xchan of one probe matches trigger chan
        vo = max(max(Q.xvalues(xp,:,bid)))./20;
        plot(X',squeeze(double(Q.xvalues(xp,:,bid)))+vo,'b');
    end
end
 text(x(end),voff,labels{2},'horizontalalignment','left','color','b','fontsize',20);
 text(x(end),0,labels{1},'horizontalalignment','left','color','k','fontsize',20);
