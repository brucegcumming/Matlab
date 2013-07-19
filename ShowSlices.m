function ShowSlices(slices, sdfs)

for j = 1:length(slices)
    [diffs, id] = min(abs(slices(j) - sdfs.times));
    sliceid(j) = id;
    sampleid(id) = id;
    delays(id) = times(id);
end
slices = sliceid;

for j = slices;
    x = [];
    y = [];
    z = [];
    sn = sampleid(j);
    for k = 1:size(sdfs.s,1)
        x(k) = sdfs.x(k,1);
        for co = 1:size(sdfs.s,2)
            if isempty(sdfs.s{k,co})
                y(k,co) = NaN;
            else
                y(k,co) = sdfs.s{k,co}(sn);
            end
        end
     end
     for co = 1:size(sdfs.s,2)
         h(np) = plot(x,y(:,co),'color',colors{np},'linestyle',linestyles{co});
         hold on;
     end
     for co = 1:size(sdfs.extras)
       plot([min(x) max(x)],[sdfs.extras{co}.sdf(j) sdfs.extras{co}.sdf(j)],'--','color',colors{np});
   end
       labels{np} = sprintf('dT = %.0fms',delays(j)/10);
       np = np+1;
end
if legendpos < 6
    legend(h,labels);
end
title(sprintf('%s',Expt.Header.Name));
end