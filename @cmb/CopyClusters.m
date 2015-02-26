function cluster = CopyClusters(cluster, ec)
%
% Copy clusters from an expt (ec) to main cluster. Only copy if set.
% that way newly defined clusters on a probe are not wiped out.

for j = 1:size(ec,1)
for k = 1:size(ec,2)
if isfield(ec{j,k},'params')
cluster{j,k} = ec{j,k};
end
end
end


