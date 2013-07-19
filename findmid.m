function mid = findmid(kernimg)

guess=size(kernimg,1)/2;
range=[(guess-5):(guess+5)];
count=0;
for i=range
    for j=range
        test=kernimg((i-5):(i+5),(j-5):(j+5));
        if ~sum(abs(diag(fliplr(test-test'))))
            count=count+1;            
            mid(count,:)=[i j];    
        end
    end
end

if count>1
    PrintMsg(0,'multiple matches');
end