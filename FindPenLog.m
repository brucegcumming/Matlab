function name = FindPenLog(monkey, pen, varargin)
%name = FindPenLog(monkey, pen, ...
%finds a penetration log file belongin to this animal and penetraion.
if nargin > 1 && isnumeric(pen)
name = ['/b/bgc/anal/' monkey '/pens/pen' num2str(pen) '.log'];
if exist(name)
    return;
else
    mycprintf('blue','%s missing\n',name);
    name= [];
end
elseif ischar(monkey)
   mnk = GetMonkeyName(monkey);
   ufl = ReadUfl(monkey);
   [a,b] = Counts(ufl.rfs(:,6));
   pen = b(1);
   name = ['/b/bgc/anal/' mnk '/pens/pen' num2str(pen) '.log'];
end
