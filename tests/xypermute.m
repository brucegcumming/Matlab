function diff =  xypermute(x,y)

xy = cat(2,x,y);

for j = 1:100
id = 1 + round(rand(size(x, 1)));
aid = cat(1, id, 2-id);

xid = sub2ind(size(xy),1:size(xy,1),id);
yid = sub2ind(size(xy),1:size(xy,1),aid);
diff(j) = mean(xy(xid)-xy(yid));
end