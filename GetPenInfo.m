function txt = GetPenInfo(T, varargin)

rid = strmatch('Right',T.text);
eid = strmatch('Electrode',T.text);
tid = strmatch('Tube',T.text);
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