function ParForIndex(loops)
%Test labindex and numlabs functions

    parfor (j = 1:loops)
        for k = 1:1000 %do enough work that it takes noticeable time
        X{j} = rand(100).^0.5;
        end
        id = labindex;
        fprintf('loop %d in Lab %d of %d\n',j,id,numlabs);
    end