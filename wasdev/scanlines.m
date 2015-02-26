function txt = scanlines(name,varargin)
%txt = scanlines(name)  retunrs a cell array of strings corresponding
%to the lines in text file name.  Wrapper for textscan+fopen
% ...,'silent', suppresses error message (Usually prints in red if file
% doen not exist)
silent = 0;
j = 1;
argon = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'bufsize',4)
        argon = {argon{:} varargin{j} varargin{j+1}};
        j = j+2;
    elseif strncmpi(varargin{j},'silent',4)
        silent = 1;
    end
    j = j+1;
end
txt = {};
if ~exist(name,'file')
    if silent == 0
        mycprintf('errors','%s Does not exist\n',name);
    end
    return;
end
fid = fopen(name,'r');
if fid < 0
    mycprintf('errors','%s exists but fopen fails\n',name);
    return;
end    
try
a = textscan(fid,'%s','delimiter','\n',argon{:});
fclose(fid);
txt = a{1};
catch
    fclose(fid);
    txt = {};
end
