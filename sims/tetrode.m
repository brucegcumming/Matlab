function [pcs, E] = tetrode(nspk, varargin)
%[pcs, E] = tetrode(nspk, varargin)
%simulate tetrode with independent Gaussian noise on each probe
%time-smoothed independtently on each probe.  returns projection onto
%pcs and the EigenVectors
%tetrode(10000000,'thr',-1.75)  shows the artifact of triggering with an OR
%rule pretty clearly
%
j = 1;
th = -0.5;
while j <= length(varargin)
    if strncmpi(varargin{j},'thre',3)
        j = j+1;
        th = varargin{j};
    end
    j = j+1;
end
V = randn([10 4 nspk]);
for j = 1:nspk
    for k = 1:4
        V(:,k,j) = smooth(V(:,k,j),5);
    end
end
ida = find(V(5,1,:) < th);
idb = find(V(5,2,:) < th);
idc = find(V(5,3,:) < th);
idd = find(V(5,4,:) < th);
id = union(ida,idb);
id = union(id,idc);
id = union(id,idd);

V = reshape(V(:,:,id), [40 length(id)]);
[E, EV] = eig(cov(V'));
pcs = V' * E;
if length(id) > 1000
    sz = 1;
else
    sz = 5;
end
subplot(2,4,1);
plot(pcs(:,end),pcs(:,end-1),'.','markersize',sz);
subplot(2,4,2);
plot(pcs(:,end),pcs(:,end-2),'.','markersize',sz);
subplot(2,4,3);
plot(pcs(:,end),pcs(:,end-3),'.','markersize',sz);
subplot(2,4,4);
plot(pcs(:,end),pcs(:,end-4),'.','markersize',sz);
subplot(2,4,5);
plot(pcs(:,end-1),pcs(:,end-2),'.','markersize',sz);
subplot(2,4,6);
plot(pcs(:,end-1),pcs(:,end-3),'.','markersize',sz);
subplot(2,4,7);
plot(pcs(:,end-1),pcs(:,end-4),'.','markersize',sz);
subplot(2,4,8);
plot(pcs(:,end-2),pcs(:,end-3),'.','markersize',sz);
