function [hostname, details] = gethostname(varargin)

os = computer;
if strncmp(os,'PCWIN',5)
    hostname=getenv('COMPUTERNAME');
    if isempty(hostname)
    hostname=getenv('USERDOMAIN');
    end
else
    [a, hostname] =system('hostname');
end
if isempty(hostname)
    hostname = 'UNKNOWNHOST';
end
details.hostid = hostid;
