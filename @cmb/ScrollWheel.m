function ScrollWheel(src, evnt)
global mousept;

DATA = GetDataFromFig(src);

if src ~= gcf
return;
end

mousept.angle = mousept.angle+0.02*evnt.VerticalScrollCount;
mousept.mode = 10;
if myhandle(mousept.lasth)
delete(mousept.lasth);
end
mousept = cmb.myellipse(mousept,[0 0; 0 0 ]);


