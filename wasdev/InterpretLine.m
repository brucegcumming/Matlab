function res = InterpretLine(txt, varargin)

% InterpretLine(s, ....)
% Reads a text string and returns extracted values according to the
% contents. By default strings are supposed to be lines send by binoc down
% the serial port.  Typically returns a number

res = [];
if isempty(txt)
    return;
end
if txt(1) == ' '
    txt = txt(2:end);
end

if sum(strncmp(txt,{'RightHem' 'Electrode'},8))
    id = strfind(txt,'Contact');
    if length(id)
        x = id(1);
        id = strfind(txt(id:end),' ');
        x = sscanf(txt(id+x:end),'%d');
        res.probesep = x;
    end
    if strncmp(txt,'Electrode',8)
        id = strfind(txt,' ');
        res.electrode = txt(1+id(1):end);
    end
    if strncmp(txt,'RightHemi',8)
        res.hemisphere = 'Right';
    end
end