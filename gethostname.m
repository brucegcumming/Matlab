function hostname = gethostname(varargin)

if ispc 
    hostname = getenv('COMPUTERNAME');
    if isempty(hostname)
        hostname=getenv('USERDOMAIN');
    end
else
    [status, hostname] = system('hostname');
end

%elseif strncmp(os,'GLNXA64',7)
%    hostname=getenv('HOSTNAME');
if isempty(hostname)
    hostname = 'UNKNOWNHOST';
end
hostname = deblank(hostname); %in case trailing cr etc