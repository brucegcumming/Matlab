function tt = TimeMark(tt, str, show)
%tt = TimeMark(tt, str, show)  keep track of exection time

if nargin < 3
    show = 0;
end
if nargin == 1  || (isnumeric(str) && str > 0)
    t = timediff([tt.time]);
    for j = 1:length(tt)
        fprintf('%s %.2f\n',tt(j).str,t(j));
    end
    return;
end

    ttn = length(tt)+1;
    tt(ttn).time = now;
    tt(ttn).str = str;
    if show
        fprintf('%s at %.2f\n',str,mytoc(tt(1).time));
    end
