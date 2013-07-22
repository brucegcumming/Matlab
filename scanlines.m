function txt = scanlines(name)
%txt = scanlines(name)  retunrs a cell array of strings corresponding
%to the lines in text file name.  Wrapper for textscan+fopen

txt = {};
if ~exist(name,'file')
    mycprintf('errors','%s Does not exist\n',name);
    return;
end
fid = fopen(name,'r');
if fid < 0
    mycprintf('errors','%s Does not exist\n',name);
    return;
end    
a = textscan(fid,'%s','delimiter','\n');
fclose(fid);
txt = a{1};
