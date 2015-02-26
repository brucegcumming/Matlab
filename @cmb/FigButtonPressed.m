function FigButtonPressed(src, data)
global mousept;
set(src, 'WindowButtonMotionFcn',@cmb.FigButtonDragged);

mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
mousept.down = 1;
mousept.lasth = 0;
mousept.drags = 0;
hold on; %othewise drawing ellipse deletes data
if mousept.mode == 1
mousept.start = get(gca,'CurrentPoint');
end
mousept.l = mousept.start;
mousept.siz = [0 0];

