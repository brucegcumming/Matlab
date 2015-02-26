function name = FindPenLog(monkey, pen, varargin)
%name = FindPenLog(monkey, pen, ...
%finds a penetration log file belongin to this animal and penetraion.

name = ['/bgc/anal/' monkey '/pens/pen' num2str(pen) '.log'];
if exist(name)
    return;
else
    mycprintf('blue','%s missing\n',name);
    name= [];
end