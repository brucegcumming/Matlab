function m = SetDiag(m, val)

for j = 1:size(m,1)
    m(j,j) = val;
end
