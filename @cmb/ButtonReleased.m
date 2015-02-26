function ButtonReleased(src, data)
global mousept;

%mousept
if mousept.down == 0 %do nothing
    if mousept.mode == 11
        mousept.mode = 0;
    end
    % delete(mousept.lasth);
    % mousept.lasth = 0;
    % Even if touching a cluster in AllProbeSpikes may need to
    % classify the spikes. Will not have happened during spooling.
    if ~strcmp(get(src,'tag'),'AllProbeSpikes')
        return;
    end
end
mousept.down = 0;
pt = get(gca,'CurrentPoint');
if myhandle(mousept.lasth)
    delete(mousept.lasth);
end
mousept= cmb.myellipse(mousept,pt)
if mousept.mode == 5 %presed inside ellise without dragging
    set(mousept.lasth,'linewidth',2);
else
    cmb.ClassifySpikes(mousept, src);
    mousept.mode = 0;
end
%mousept.drags

