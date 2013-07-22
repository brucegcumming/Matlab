function txt = FixMat(txt, res)

frid = strmatch('Fr',txt.text);
for j = 1:length(res.id)
    idstr = sprintf('id%d',res.id(j));
    id = strmatch(idstr,txt.text);
    if length(id) == 1
        if strncmp('cs',txt.text(id-1,:),2)
            str = sprintf('Fr%.2f',res.Fr(j));
            txt.text(id-1,1:length(str)) = str;
        end
    end
end