function DATA = LoadComments(DATA, dirname)

outname = [dirname '/Comments.mat'];
if exist(outname,'file')
    load(outname);
    if exist('Comments','var')
        DATA.Comments = Comments;
    end
    if exist('Tagged','var') && ~isempty(Tagged)
        DATA.tagged = Tagged;
    end
end