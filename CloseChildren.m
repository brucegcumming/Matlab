function CloseChildren(parent, varargin)
%CloseChildren(parent, closes figures associated with parent
% - idendifed via appdata'ParentFigure' in the Childrens

closeall = 0;
j = 1;
args = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'closeall',8)
        closeall = 1;
    end
    j = j+1;
end
%find hidden figure too
%figs = findobj('type','figure');
c = allchild(0);
types = get(c,'type');
fid  = find(strcmp('figure',types));
figs = c(fid);
figs = setdiff(figs,parent);
for j = 1:length(figs)
    if getappdata(figs(j),'ParentFigure') == parent
        if ~isappdata(figs(j),'KeepFigure') || closeall
            if figs(j) ~= parent
                close(figs(j));
            end
        end
    end
end
