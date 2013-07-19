function SaveComments(DATA, dirname)

outname = [dirname '/Comments.mat'];
Comments = DATA.Comments;
if isfield(DATA,'tagged')
    Tagged = DATA.tagged;
else
    Tagged = [];
end
save(outname,'Tagged','Comments');
