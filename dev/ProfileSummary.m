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

for j = 1:length(P.FunctionTable)
    F = P.FunctionTable(j);
    for k =1:length(fclist)
        if strfind(F.FunctionName,fclist{k})
            ListLines(P.FunctionTable, j);
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
    fprintf('Function %s (%d)  %.3f(%.3f):\n',F.FunctionName,id,a(j),sum(F.ExecutedLines(:,3)));
end
fprintf('Total selftime %.3f\n',sum(selftime));


function ListLines(T, id)

F = T(id);

fprintf('Function %s (%d)  %.3f(%.3f):\n',F.FunctionName,id,F.TotalTime-sum([F.Children.TotalTime]),sum(F.ExecutedLines(:,3)));
[a,b] = sort(F.ExecutedLines(:,3));
b = b(a> 0);
for j = 1:length(b)
    fprintf('%s %d %d %.3f\n',F.FunctionName,F.ExecutedLines(b(j),1),F.ExecutedLines(b(j),2),F.ExecutedLines(b(j),3));
end
fprintf('Children:\n');
for j = 1:length(F.Children)
    fprintf('%s %.3f\n',T(F.Children(j).Index).FunctionName,F.Children(j).TotalTime);
end
fprintf('Parents:\n');
for j = 1:length(F.Parents)
    fprintf('%s %d\n',T(F.Parents(j).Index).FunctionName,F.Parents(j).NumCalls);
end



