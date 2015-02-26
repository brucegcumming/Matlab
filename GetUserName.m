function username = GetUserName(varargin)

verbose = 0;
j = 1;
while j <= length(varargin)
    if sum(strncmpi(varargin{j},{'verbose' 'show' 'print'},4))
        verbose = 1;
    elseif strncmpi(varargin{j},'force',4)
        j = j+1;
        username = varargin{j};
        setappdata(0,'UserNameForce',username);
        return;
    end
    j = j+1;
end

username = getappdata(0,'UserNameForce');
if ~isempty(username)
    username = [username '*'];
    return;
end
    
os = computer;
rootpath = 'Z:';
if strncmp(os,'PCWIN',5)
    username=getenv('username');
    rootpath = 'Z:';
elseif strncmp(os,'GLNXA64',7)
    username=getenv('USER');
elseif strncmp(os,'MACI64',6)
    username=getenv('USER');
    rootpath = '/b';
else
    username=getenv('USER');
    rootpath = '/b';
end


if sum(strcmpi(username,{'Expt' 'bgcgroup'}))
    usernames = {'bgc' 'bondya' 'moeenya' 'kangi' 'sid' 'aparicio'};
    setappdata(0,'ExptUserName',username);
    for j = 1:length(usernames)
        fid = fopen(sprintf('%s/%s/testmatlabwrite',rootpath,usernames{j}),'w');
        if fid > 0
            if verbose
                fprintf('Hello Expt, I think you are %s (from %s)\n',usernames{j},mfilename('fullpath'));
            end
            username = [username ' (' usernames{j} ')'];
            setappdata(0,'ExptUserName',usernames{j});
            fclose(fid);
        break;
        end
    end
end