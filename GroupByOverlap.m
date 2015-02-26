function [igroups, details] = GroupByOverLap(X, varargin)
%Given a matix of mutual overlaps, break into groups


    igroups = {};
    ogroups = {};
    groups = {};
%Anywhere there is nothing on the diagonal 2 away from the center, it
%sugests a different grouping
     diaga(1) = 1;
     len = size(X,1);
     for j = 2:size(X,1)
         diaga(j) = min([X(j,j-1) X(j-1,j)]);
     end
     for j = 1:(2*size(X,1) -3)
         o = ceil(j/2);
         p = 2 + floor(j/2);
         ndiag(j,1) = X(p,o);
         ndiag(j,2) = X(o,p);
         ppos(j,1) = o;
         ppos(j,2) = p;         
     end
     breakcrit = mean(ndiag');
     inid = find(breakcrit > 0.5);
     outid = find(breakcrit < 0.1);
     if isempty(inid)
         ingroup = 0;
         if isempty(outid)
             outid = min(breakcrit);
         end
         j = outid(1);
         group = [];
     elseif isempty(outid) || inid(1) < outid(1)
         ingroup = 1;
         group = [inid(1) inid(1)+1];
         j = inid(1);
     else
         ingroup = 0;
         j = outid(1);
         group = [];
     end
         pregroup = [];
         ogroup = [];
     ng=0;
     while j <= length(breakcrit)
         if isempty(group) 
             if breakcrit(j) > 0.5
                 group = ppos(j,:);
                 pregroup = setdiff(unique(pregroup),group);
                 for p = pregroup(:)'
                     x = [X(p,group); X(group,p)'];
                     if min(mean(x')) > 0.1
                         group = [group p];
                     else
                         ogroup = [ogroup p];
                     end
                 end
             else
                 pregroup = [pregroup ppos(j,:)];
             end
         else
             if breakcrit(j) < 0.1
                 ng = ng+1;
                 igroups{ng} = group;
                 ogroups{ng} = ogroup;
                 group = [];
                 pregroup = [];
             elseif breakcrit(j) > 0.5
                 group = unique([group ppos(j,:)]);
             else
                 p = setdiff(ppos(j,:),group);
                 if ~isempty(p)
                 x = [X(p,group); X(group,p)'];
                 if min(mean(x')) > 0.1
                     group = [group p];
                 else
                     ogroup = [ogroup p];
                 end
                 end
             end
         end
         j = j+1;
     end
     if ~isempty(group)
         ng = ng+1;
         igroups{ng} = group;
         ogroups{ng} = ogroup;
     end
     x = cat(2,igroups{:});
     for j = 1:length(ogroups)
         ogroups{j} = unique(ogroups{j});
         if ~isempty(ogroups{j})
         x = cat(2,x,ogroups{j});
         end
     end
     details.nogroup = setdiff(1:size(X,1),x);
     details.ogroups = ogroups;
