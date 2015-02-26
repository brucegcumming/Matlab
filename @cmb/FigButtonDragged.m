function FigButtonDragged(src, data)
global mousept;

if mousept.down
pt = get(gca,'CurrentPoint');
if mousept.lasth & ishandle(mousept.lasth)
delete(mousept.lasth);
end
mousept= cmb.myrect(mousept,pt);
%   plot(pt(1,1),pt(2,2),'+');
mousept.drags = mousept.drags+1;
end

