function txt = scanlines(name,varargin)
%txt = scanlines(name)  retunrs a cell array of strings corresponding
%to the lines in text file name.  Wrapper for textscan+fopen
% ...,'silent', suppresses error message (Usually prints in red if file
% doen not exist)
%scanlines (...,'matrix') scans each lines into doubles, and builds a
%matrix. Does not stop if a line is short, but Nan pads;
silent = 0;
mkmatrix = 0;
method = 0;
j = 1;
argon = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'bufsize',4)
        argon = {argon{:} varargin{j} varargin{j+1}};
        j = j+2;
    elseif strncmpi(varargin{j},'matrix',4)
        mkmatrix = 1;
    elseif strncmpi(varargin{j},'fget',4) %use fgets instead of textread
        method = 1;
    elseif strncmpi(varargin{j},'silent',4)
        silent = 1;
    elseif ischar(varargin{j})
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

if method == 1
    txt = {};
    line = fgets(fid);
    while line > 0
        txt{end+1} = deblank(line); %remove \n
        line = fgets(fid);
    end
    fclose(fid);
else
    
try
a = textscan(fid,'%s','delimiter','\n',argon{:});
fclose(fid);
txt = a{1};
catch
    fclose(fid);
    txt = {};
end
end

if mkmatrix
    m = [];
    for j = 1:length(txt)
        f = sscanf(txt{j},'%f');
        if ~isempty(f)
            m(end+1,1:length(f)) = f;           
            if length(f) < size(m,2)
                m(end,length(f)+1:end) = NaN;
            end
        end
    end
    txt = m;
end
