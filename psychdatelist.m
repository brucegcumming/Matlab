function datelist = psychdatelist()
k=0;
for a=2008:2011
    for b=1:12
        for c=1:32
            k=k+1;
            datelist(k,:)=[b-1,c-1,a];
        end
    end
end