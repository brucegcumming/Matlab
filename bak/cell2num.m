function x = cell2num(C,f)

for j = 1:length(C)
    for k = 1:length(f)
        x(j,k) = eval(['C{j}.' f{k}]);
    end
end