function r = hostid()
if ispc 
    r = getenv('COMPUTERNAME');
else
    r = system('hostname');
end