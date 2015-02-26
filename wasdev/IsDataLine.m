function yn = IsDataLine(h, varargin)

if ~isempty(h) && ishandle(h(1)) && strcmp(get(h(1)),'line');
    yn = 1;
else
    yn = 0;
end