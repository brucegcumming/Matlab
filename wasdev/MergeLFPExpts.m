function [xcs, bxc, lfpall, rc] = MergeLFPExpts(LFP,varargin)

%MergeLFPExpts(LFP)
%Combine expts run at different depths into one large
%map. Calculates a map for each block, and finds vertical displacement that
%maximizes correlation.
blocks = 1:length(LFP.Header.BlockStart);
offsets = zeros(size(blocks));
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'blocks',5)
        j = j+1;
        blocks = varargin{j};
    elseif strncmpi(varargin{j},'offsets',5)
        j = j+1;
        offsets = varargin{j};
    end
    j = j+1;
end

for j = 1:length(blocks)
    blk = blocks(j);
    rc = PlotRevCorAny(LFP,'lfp','block',blk,'nmin',10);
    lfpall(:,:,:,:,j) = RC2mat(rc);
    lfpn(j) = mean(rc.sdfs.n(:));
    xv(:,j) = rc.sdfs.x(:,1);
    if ~isempty(rc.sdfs.extras)
    blanks(:,:,j) = rc.sdfs.extras{1}.lfp;
    end
end

subplot(2,size(lfpall,5),1);
imagesc(squeeze(var(lfpall(:,:,1,:,1),[],2))');
cx(1,:) = caxis;
title(sprintf('%.2f',LFP.Header.depths(1)));
nch = size(lfpall,4);
nex = size(lfpall,5);
mresp = mean(mean(lfpall,2),3);
nx = size(lfpall,2);
ny = size(lfpall,3);
mresp = repmat(mresp,[1 nx ny 1 1]);
for j = 1:size(blanks,3)
    blanks(:,:,j) = blanks(:,:,j) - squeeze(mresp(:,1,1,:,j));
end
lfpall = lfpall - mresp;
for blk = 2:size(lfpall,5)
subplot(2,size(lfpall,5),blk);
imagesc(squeeze(var(lfpall(:,:,1,:,blk),[],2))');
title(sprintf('%.2f',LFP.Header.depths(blk)));
cx(blk,:) = caxis;
for j = -22:0
    a = lfpall(:,:,:,1:nch+j,blk-1);
    b = lfpall(:,:,:,1-j:nch,blk);
    xc = corrcoef(a(:),b(:));
    xcs(j+23,blk) = xc(1,2);
    diffs(j+23,blk) = mean(abs(a(:)-b(:)));
end
for j = 1:22
    a = lfpall(:,:,:,1+j:nch,blk-1);
    b = lfpall(:,:,:,1:nch-j,blk);
    xc = corrcoef(a(:),b(:));
    xcs(j+23,blk) = xc(1,2); 
    diffs(j+23,blk) = mean(abs(a(:)-b(:)));
end
end

for blk = 1:size(lfpall,5)
subplot(2,size(lfpall,5),blk);
caxis([min(cx(:)) max(cx(:))]);
end


subplot(2,nex,nex+1);
imagesc(squeeze(blanks(:,:,1))');
bx(:,1) = caxis;
for blk = 2:size(blanks,3)
subplot(2,nex,nex+blk);
imagesc(squeeze(blanks(:,:,blk))');
bx(:,blk) = caxis;
for j = -22:0
    a = blanks(:,1:nch+j,blk-1);
    b = blanks(:,1-j:nch,blk);
    xc = corrcoef(a(:),b(:));
    bxc(j+23,blk) = xc(1,2);
    bdiffs(j+23,blk) = mean(abs(a(:)-b(:)));
end
for j = 1:22
    a = blanks(:,1+j:nch,blk-1);
    b = blanks(:,1:nch-j,blk);
    xc = corrcoef(a(:),b(:));
    bxc(j+23,blk) = xc(1,2);
    bdiffs(j+23,blk) = mean(abs(a(:)-b(:)));
end
end

for blk = 1:size(blanks,3)
subplot(2,nex,blk+nex);
caxis([min(bx(:)) max(bx(:))]);
end

sz = size(lfpall);
sz(4) = sz(4) + max(offsets);
merged = zeros(sz(1:4));
nx = size(lfpall,5);
mn = zeros(nx,sz(4));
for blk = 1:size(lfpall,5)
    merged(:,:,:,[1:23]+offsets(blk)) = merged(:,:,:,[1:23]+offsets(blk)) + lfpall(:,:,:,1:23,blk).*lfpn(blk);
    mn(blk,[1:23]+offsets(blk)) = mn(blk,[1:23]+offsets(blk)) + lfpn(blk); 
end
mn = sum(mn);
for j = 1:size(merged,4)
    merged(:,:,:,j) = merged(:,:,:,j)./mn(j);
end

sz = size(blanks);
sz(2) = sz(2)+max(offsets);
bmerged = zeros(sz(1:2));
mn = zeros(nx,sz(2));
for blk = 1:size(blanks,3)
    bmerged(:,[1:23]+offsets(blk)) = bmerged(:,[1:23]+offsets(blk)) + blanks(:,1:23,blk).*lfpn(blk);
    mn(blk,[1:23]+offsets(blk)) = mn(blk,[1:23]+offsets(blk)) + lfpn(blk); 
end
mn = sum(mn);
for j = 1:size(bmerged,4)
    bmerged(:,:,:,j) = bmerged(:,:,:,j)./mn(j);
end

rc.lfp = merged;
rc.lfpblank =bmerged;
rc.x = rc.sdfs.x(:,1);
rc.y = rc.sdfs.y(1,:);
rc.lfptimes = rc.sdfs.lfptimes;
rc.isrc = 1;
