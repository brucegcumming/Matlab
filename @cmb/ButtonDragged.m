function ButtonDragged(src, data)

global mousept;

if ~isfield(mousept,'down')
    return;
end
%fprintf('mode %d %.4f,%d\n',mousept.mode,mousept.angle,mousept.drags);
if mousept.down && mousept.mode == 5
    mousept.drags = mousept.drags+1;
    if mousept.drags > 2
        mousept.down = 1;
        mousept.mode = 11;
    end
elseif mousept.down
    pt = get(gca,'CurrentPoint');
    if myhandle(mousept.lasth)
        delete(mousept.lasth);
    end
    mousept= cmb.myellipse(mousept,pt);
    %    plot(pt(1,1),pt(2,2),'+');
    mousept.drags = mousept.drags+1;
end


