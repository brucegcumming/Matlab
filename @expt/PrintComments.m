function PrintComments(Expt,varargin)
Comments = [];
if iscell(Expt)
    for j = 1:length(Expt)
       C{j} = expt.PrintComments(Expt{j});
    end
    Comments = cat(1,C{:});
    return;
end

if ~isfield(Expt,'Comments')
    fprintf('No Comments in %s\n',GetName(Expt));
    return;
end
for j = 1:length(Expt.Comments)
    t = ConvertTime(Expt, Expt.Comments.times(j));
    fprintf('%s:%s\n',datestr(t),Expt.Comments.text{j});
    Comments(j).text = Expt.Comments.text{j};
    Comments(j).times = t;
end