function good = GoodCluster(C, varargin)

crit(1,:) = [1 0.015 0.25 2];
crit(2,:) = [2 0.015 0.4 3];
good = 0;
if ~isfield(C,'dipsize')
    C.dipsize = 0;
end
    
if C.dipsize(1) > crit(2,1) & C.hdip > crit(2,1) || C.bmc > crit(2,3) || C.mahal(1) > crit(2,4);
    good = 2;
elseif C.dipsize(1) > crit(1,1) || C.hdip > crit(1,2) || C.bmc > crit(1,3) || C.mahal(1) > crit(1,4)
    good = 1;
end