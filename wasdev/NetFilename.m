function [netname, netdir] = NetFilename(name, varargin)

[a,b] = fileparts(name);
monkey = GetMonkeyName(name);
name = strrep(name,'\','/')
netname = regexprep(name, ['.*/' monkey '/'],['Z:/bgc/data/' monkey '/']);
if strcmp(netname, name) %no match
    netname = [];
    netdir = 0;
else
    [a,b] = fileparts(netname);
    netdir = a;
end