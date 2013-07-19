function binoc(varargin)


go = 1;
win = figure('Tag','BinocLoop');

gobtn = uicontrol(gcf,'Style', 'CheckBox','String','Go','Units', 'normalized','Position', [0.1 0.1 0.3 0.3],...
      'Tag','GoButton','value',1);

drawnow;

warning('off','MATLAB:serial:fscanf:unsuccessfulRead');
port = serial('COM1', 'BaudRate', 1200)
set(port,'Databits',7);
set(port,'Parity','even');
set(port,'StopBits',1);
set(port,'FlowControl','hardware');
set(port,'Terminator','CR/LF');
set(port,'Timeout',0.01);
fopen(port);
while go
    line = fscanf(port);
    if length(line)
        fprintf('%s\n',line);
    end
    go = get(gobtn,'value');
end
close(port);