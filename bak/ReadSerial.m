function ReadSerial(varargin)

topw = makewindow;
drawnow;
pause(0.1);

portname = 'COM1';
baudrate = 115200;

nlines = 10;
j = 1;
while j < nargin+1
    if strncmpi(varargin{j},'baud',4)
        j = j+1;
        baudrate = varargin{j};
    end
    j = j+1;
end

fprintf('Using %d baud\n',baudrate);
s = serial(portname, 'BaudRate', baudrate);
set(s,'Databits',7);
set(s,'Parity','even');
set(s,'StopBits',1);
set(s,'FlowControl','hardware');
set(s,'Terminator','CR/LF');

stop = 0;

fopen(s);
j = 0;
while j < nlines & ~stop
    line = fscanf(s)
    if ~isempty(line)
        fprintf(s,'Got that\n');
    end
    stop = get(findobj('Tag','StopSerial'),'value');
    j = j+1;
end
fclose(s);
stop = 0;
close(topw);

function cntrl_box = makewindow(varargin)

name = 'ReadSerial';
tag = name;

    
cntrl_box = figure('Position', [10 10 100 40],...
    'NumberTitle', 'off', 'Tag',tag,'Name',name);

uicontrol(gcf,'Style', 'checkbox','String', 'Stop', 'Tag', 'StopSerial', 'Position', [10 10 50 20]);