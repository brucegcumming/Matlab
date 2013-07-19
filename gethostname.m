function [hostname, details] = gethostname(varargin)

os = computer;
if strncmp(os,'PCWIN',5)
    hostname=getenv('COMPUTERNAME');
    if isempty(hostname)
    hostname=getenv('USERDOMAIN');
    end
elseif strncmp(os,'GLNXA64',7)
    hostname=getenv('HOSTNAME');
else
    [a, hostname] =system('hostname');
end
if isempty(hostname)
    hostname = 'UNKNOWNHOST';
end
details.hostid = hostid;
hostname = deblank(hostname); %in case trailing cr etc