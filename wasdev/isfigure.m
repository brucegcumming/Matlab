function yesno = isfigure(h)
% function yesno = isfigure(h)% Returns 1 if the handle h is a figure, and 0 otherwise% % Jenny Read, 03/14/2003
if ishandle(h)
    yesno = ~isempty(findobj(h,'flat','type','figure'));else    yesno=0;end