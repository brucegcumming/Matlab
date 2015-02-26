function UpdateCode(varargin)
listonly = 0;
args = {};
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'list',4)
    listonly = 1;
end
j = j+1;
end
if listonly == 0
    args = {args{:} 'copynewer'};
end
        
    
cpfiles('/b/bgc/matlab','/b/bgc/matlab/AllV',args{:},'hidenew');
