function val = getstr2num(tag)

% val = getstr2num(tag) 
% reads the contents of a text field in a gui (with Tag tag)
% and returns the contents as a number
%

val = NaN;
te = findobj('Tag',tag);
if ~isempty(te)
    str = get(te(1),'String');
    val = str2num(str);
end

if length(te) > 1
    warning('Too many tags');
end