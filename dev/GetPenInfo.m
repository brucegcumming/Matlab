function txt = GetPenInfo(T, varargin)

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'name',4)
    end
    j = j+1;
end

if ~iscellstr(T.text) %worth converting for faster GetPenInfo
    T.text = CellStr(T.text);
end
rid = find(strncmp('Right',T.text,5));
eid = find(strncmp('Electrode',T.text,9));
if isempty(eid)
    eid = find(strncmp(' Electrode',T.text,9));
end
tid = find(strncmp('Tube',T.text,4);
allid = sort([rid; eid; tid]);
if isempty(allid)
    id = strmatch('uf',T.text(1:10,:));
    if length(id)
        name = regexprep(T.text(id(1),:),'.*\\','');
        name = regexprep(name,'\n.*','');
        txt = ['No Pen info in ' deblank(name)];
    else
        txt = 'No Pen info in File';
    end
    return;
end


if length(eid) && length(rid)
    txt = sprintf('%s %s',deblank(T.text(rid(end),:)),deblank(T.text(eid(end),:)));
else 
    txt = deblank(T.text(allid(end),:));
end