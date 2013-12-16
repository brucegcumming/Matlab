function ProfileSummary(P, varargin)
% ProfileSummary(P, varargin)  summarizes results from P = profile('info')
%mostly useful for  version where the profile viewer is broken (!! eg 2012b)

fclist = {};

sortby = 'self';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose',6)
    elseif strncmpi(varargin{j},'elapsed',6)
        sortby = 'total';
    elseif ischar(varargin{j})
        fclist = {fclist{:} varargin{j}};
    end
    j = j+1;
end

filenames = unique({P.FunctionTable.FileName});
P.src{length(filenames)} = [];
for j = 1:length(P.FunctionTable)
    P.FunctionTable(j).fileid = find(strcmp(P.FunctionTable(j).FileName,filenames));
end
for j = 1:length(P.FunctionTable)
    F = P.FunctionTable(j);
    for k =1:length(fclist)
        if strfind(F.FunctionName,fclist{k})
            P = LoadSrc(P, j);
            ListLines(P, j);
        end
    end
    selftime(j) = F.TotalTime - sum([F.Children.TotalTime]);
end
if strcmp(sortby,'total')
    [a,b] = sort([P.FunctionTable.TotalTime],'descend');
else
    [a,b] = sort(selftime,'descend');
end
for j = 1:min([10 length(b)])
    id = b(j);
    F = P.FunctionTable(b(j));
    fprintf('Function %s (id %d, %d calls)  %.3f(%.3f):\n',F.FunctionName,id,F.NumCalls,a(j),sum(F.ExecutedLines(:,3)));
end
fprintf('Total selftime %.3f\n',sum(selftime));

function P = LoadSrc(P, id)
   srcid = P.FunctionTable(id).fileid;
   if isempty(P.src{srcid})
       P.src{srcid} = scanlines(P.FunctionTable(id).FileName);
   end

function ListLines(P, id)

T = P.FunctionTable;
F = T(id);
src = P.src{F.fileid};

fprintf('Function %s (%d)  %.3f(+%.3f=%.3f):\n',F.FunctionName,id,F.TotalTime-sum([F.Children.TotalTime]),sum([F.Children.TotalTime]),sum(F.ExecutedLines(:,3)));
[a,b] = sort(F.ExecutedLines(:,3));
b = b(a> 0);
for j = 1:length(b)
    fprintf('%s %d %d %.3f',F.FunctionName,F.ExecutedLines(b(j),1),F.ExecutedLines(b(j),2),F.ExecutedLines(b(j),3));
    if ~isempty(src)
        fprintf('\t%s',src{F.ExecutedLines(b(j),1)});
    end
    fprintf('\n');
end
fprintf('Children:\n');
for j = 1:length(F.Children)
    fprintf('%s %.3f\n',T(F.Children(j).Index).FunctionName,F.Children(j).TotalTime);
end
fprintf('Parents:\n');
for j = 1:length(F.Parents)
    fprintf('%s %d\n',T(F.Parents(j).Index).FunctionName,F.Parents(j).NumCalls);
end



