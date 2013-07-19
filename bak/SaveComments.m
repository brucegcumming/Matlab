function SaveComments(DATA, dirname)

outname = [dirname '/Comments.mat'];
Comments = DATA.Comments;
Tagged = DATA.tagged;
save(outname,'Tagged','Comments');
