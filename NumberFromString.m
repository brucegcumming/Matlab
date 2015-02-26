function val = NumberFromString(str, pat)
%find pattern in string, and get following number

id = strfind(str,pat);
if isempty(id)
    val = [];
else
    pos = id(1)+length(pat);
    val = sscanf(str(pos:end),'%f');
    if isempty(val) && ismember(str(pos),'= ')
        pos = pos+1;
        val = sscanf(str(pos:end),'%f');
    end
end