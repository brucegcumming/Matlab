
x = instrfind('type','serial');
delete(x);
s = serial('COM1','BaudRate',9600);
fopen(s);
[F, isnew] = GetFigure('SerialPort');
if isnew
    stoph = uicontrol('style','check','string','Stop','tag','StopToggle');
    drawnow;
else
    stoph = findobj(F,'tag','StopToggle');
    set(stoph,'value',0);
end
stop = 0;
while stop == 0
    sin = fscanf(s);
    if length(sin)
        fprintf('%s\n',sin);
        if strfind(sin,'LA');
            id = strfind(sin,'LA');
            newpos = sscanf(sin(id+2:end),'%d')
        elseif strfind(sin,'POS');
            fprintf('100');
            fprintf(s,'100\n')
        end
    end
    stop = get(stoph,'value');
end