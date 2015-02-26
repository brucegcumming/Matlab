function name = BuildName(name)

if isempty([strfind(name,'/') strfind(name,'\')])
name = [pwd '/' name];
end
id = strfind(name,':');
if id
if isunix
name = ['/bgc' name(id(1)+1:end)];
name = strrep(name,'\','/');
else
name = name(id(1)+1:end);
end
else
name = name;
end


