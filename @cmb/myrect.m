function mousept = myrect(mousept, finish)
start = mousept.start;
mousept.finish = finish;
mousept.siz = finish-start;
if mousept.mode == 1 % horizontal line
mousept.lasth = plot([start(1,1) finish(1,1)],[start(2,2) start(2,2)],'r');
else
mousept.lasth = plot([start(1,1) finish(1,1)],[start(2,2) finish(2,2)],'r');
end

