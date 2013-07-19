function TestFullV(FullV, varargin);

ispk = 1;
th = -0.5;
spts = -8:32;
V= FullV.V;
tic;
sgn = diff(sign(diff(V(ispk,:),1,2)),1,2);
 id = find(sgn > 0 & V(ispk,2:end-1) < th);
 toc;
 tic;
allid = repmat(id,length(spts),1) + repmat(spts',1,length(id));
AllV = reshape(V(ispk(1),allid),[size(allid)]);
toc
V = FullV.V';
tic;
sgn = diff(sign(diff(V(:,ispk))));
 id = find(sgn > 0 & V(2:end-1,ispk) < th);
 toc; tic
allid = repmat(id,1,length(spts)) + repmat(spts,length(id),1);
AllV = reshape(V(allid,ispk(1)),[size(allid)]);
toc
