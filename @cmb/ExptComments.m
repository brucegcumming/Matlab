function s = ExptComments(Expt, varargin)
verbose = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'all',3)
        verbose = 2;
    elseif strncmpi(varargin{j},'verbose',6)
        verbose = 1;
    end
    j = j+1;
end
if isfield(Expt,'Comments') & ~isempty(Expt.Comments.text);
    s = Expt.Comments.text;
    bid = find(strncmp('cm=back=',s,8) | strncmp('cm=noback',s,8));
    if verbose == 0
        bid = [bid find(strncmp('Experimenter',s,8) | strncmp('Electrode',s,8))];
    end
    gid = setdiff(1:length(s),bid);
    for j = gid
        fprintf('%s\n',s{j});
    end
else
    s = {};
end

