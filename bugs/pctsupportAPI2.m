function pctsupportAPI2(varargin)
allProfiles = parallel.clusterProfiles();
defaultProfile = parallel.defaultClusterProfile();
x = onCleanup(@() diary('off'));
if nargin == 0
    % select the default configuration
    config = defaultProfile;
    port = -1;
elseif nargin == 1
    config = varargin{1};
    port = -1;
elseif nargin == 2
    config = varargin{1};
    port = varargin{2};
end

if exist('pctsupport.txt','file')
    delete('pctsupport.txt')
end
diary pctsupport.txt

% Gather local system information
disp('===============Client Information==================')
iGatherLocalSystemInformation();
% Display all configurations on the system
allProfiles = parallel.clusterProfiles();
defaultProfile = parallel.defaultClusterProfile();
disp('There are following cluster profiles on this system')
fprintf('\t%s\n\n', allProfiles{:});
fprintf('%s\n\n',['The default cluster profile is: ', defaultProfile]);
disp('============Beginning Cluster Profile Tests===========')
fprintf('%\n\n',['Testing ',config,' profile']);

% Test scheduler
success = iTestParCluster(config);
if ~success
    disp('Failed to create a cluster object')
    return
end

% Test distributed jobs
success = iTestIndependentJobs(config);
if ~success
    disp('Failed running a distributed job');
    return
end
% Gather information about the worker
success = iGatherWorkerConfiguration(config);
if ~success
    disp('Failed gathering worker configuration');
    return
end

% Test connectivity back to the client
success = iTestConnectivityBackToClient(config, port);
if ~success
    disp('Connectivity Tests failed')
    return
end
%{
    %% Test parallel jobs
    success = iTestParallelJobs(config);
    if ~success
        disp('Unable to run parallel jobs');
        continue
    end
%% Test matlabpool
    success = iTestMatlabpoolJobs(config);
    if ~success
        disp('Failed MATLABPOOL tests')
        continue
    end
%}

diary off;
end

% gatherLocalHostNames
% Return local hostnames and port ranges as reported by pctconfig,
% JAVA and system('hostname')
function [hostnames, portrange] = gatherLocalHostNames()
% Gather hostname and portrange as reported by pctconfig
output = pctconfig();
hostnames{1} = output.hostname;
portrange = output.portrange;
disp('MATLAB determined the following hostname and port range');
disp(['Hostname: ', output.hostname]);
disp(['Portrange: ', num2str(output.portrange)]);

% Gather hostname information as reported by JAVA
disp('JAVA has determined the following hostnames for the machine');
jHostname = java.net.InetAddress.getLocalHost.getHostName();
hostnames{2} = char(jHostname);
disp(['Hostname: ' char(jHostname)])

jHostname = java.net.InetAddress.getLocalHost.getCanonicalHostName();
hostnames{3} = char(jHostname);
disp(['Full Hostname: ' char(jHostname)])

jHostname = java.net.InetAddress.getLocalHost.getHostAddress();
hostnames{4} = char(jHostname);
disp(['IP Address: ' char(jHostname)]);

% Gather hostname information as reported by system('hostname')
[a, b] = system('hostname');
hostnames{5} = strtrim(b);
disp(['Hostname as returned by: system(''hostname'') ', b]);
end

% iGatherWorkerConfiguration
% Submit a distributed job to gather worker configuration information
function ok = iGatherWorkerConfiguration(config)
timeOutValue = 600; %seconds
myCluster = parcluster(config);
job = createJob(myCluster);
task = createTask(job, @iGatherLocalSystemInformation, 0, {});
task.CaptureDiary = true;
job.AttachedFiles = {'pctsupportAPI2.m'};
submit(job)
% test pct connectivity
% wait for the job to finish
ok = wait(job,'finished',timeOutValue);
if ~ok
    disp('We timed out while waiting for job to run on the cluster')
    disp('Consider increasing the time out value')
    % We timed out destroy the job
    job.delete
    return
end
disp('===============Worker Information==================')
disp('----------------------------------------------------------')
disp(job.Tasks(1).Diary)
disp('----------------------------------------------------------')

end

% iTestconnectivityBackToClient
% Submit a distributed job and attempt to connect back to the client given
% the hostnames and ip address from pctconfig, JAVA and system('hostname')
function success = iTestConnectivityBackToClient(config, port)
[hostnames, portrange] = gatherLocalHostNames();
if port == -1
    port = portrange(1) + ceil((portrange(2) - portrange(1))/2);
end
try
    success = true;
    myCluster = parcluster(config);
    
    % Test with hostname returned from pctconfig
    disp('Testing connectivity with hostname returned by pctconfig')
    disp('----------------------------------------------------------')
    ok = pItestConnectivityBackToClient(myCluster, hostnames{1}, port);
    if ok
        fprintf('%s\n\n','Successfully connected back to client using hostname provided by pctconfig')
    else
        fprintf('%s\n\n','Failed to connect back to client using hostname provided by pctconfig')
        success = false;
    end
    
    % Test with hostname returned by java
    
    disp('Testing connectivity with hostname returned by java')
    disp('----------------------------------------------------------')
    ok = pItestConnectivityBackToClient(myCluster, hostnames{2}, port);
    if ok
        fprintf('%s\n\n','Successfully connected back to client using hostname provided by java')
    else
        fprintf('%s\n\n','Failed to connect back to client using hostname provided by java')
        success = false;
    end
    
    % Test with full hostname returned by java
    disp('Testing connectivity with full hostname returned by java')
    disp('----------------------------------------------------------')
    ok = pItestConnectivityBackToClient(myCluster, hostnames{3}, port);
    if ok
        fprintf('%s\n\n','Successfully connected back to client using full hostname provided by java')
    else
        fprintf('%s\n\n','Failed to connect back to client using full hostname provided by java')
        success = false;
    end
    
    % Test with IP address returned by java
    disp('Testing connectivity with IP address returned by java')
    disp('----------------------------------------------------------')
    ok = pItestConnectivityBackToClient(myCluster, hostnames{4}, port);
    if ok
        fprintf('%s\n\n','Successfully connected back to client using full hostname provided by java')
    else
        success = false;
        fprintf('%s\n\n','Failed to connect back to client using full hostname provided by java')
    end
    
    % Test with hostname returned by system('hostname')
    disp('Testing connectivity with hostname returned by system(''hostname'')')
    disp('----------------------------------------------------------')
    ok = pItestConnectivityBackToClient(myCluster, hostnames{5}, port);
    if ok
        fprintf('%s\n\n','Successfully connected back to client using hostname returned by system(''hostname'')')
    else
        success = false;
        fprintf('%s\n\n','Failed to connect back to client using hostname returned by system(''hostname'')')
    end
    disp('----------------------------------------------------------')
    
catch err
    err.getReport('extended')
    success = false;
end
end

% pItestConnectivityBackToClient
% Creates a distributed job on specified scheduler and specifies the
% hostname and port to use by the worker to attempt to connect back to
% client
function status = pItestConnectivityBackToClient(myCluster, hostname, port)
timeOut = 240; %seconds
disp(['Attempting to test connectivity with hostname = ', hostname]);
job = createJob(myCluster);
task = createTask(job, @remoteTestConnectivityBackToClinet, 1, {hostname, port});
task.CaptureDiary = true;
job.AttachedFiles = {'pctsupportAPI2.m'};
submit(job)
% test whether worker can establish connection back to client
status = establishServerSocket(port);
if(status)
    disp('successfully received connection from a worker');
else
    disp('failed to receive a connection from worker');
end

% wait for the job to finish
completed = wait(job,'finished',timeOut);
if ~completed
    % we timed out destroy the job
    disp('Job Timed Out')
    disp('----------------------------------------------------------')
    job.Tasks(1)
    disp(job.Tasks(1).Diary)
    job.delete;
else
    disp('Information about the task')
    disp('----------------------------------------------------------')
    job.Tasks(1)
    fprintf('Task Command Window Output: \n\n')
    disp(job.Tasks(1).Diary)
end
end

% remoteTestConnectivityBackToClient
% This function executes on the worker and attempts to connect back the
% client using java sockets given a specific hostname and port.
% Hostname can actually be an IP address. This function also tries to
% perform nslookup and ping using the client hostname from the worker.
function ok = remoteTestConnectivityBackToClinet(hostname, port)
[status, result] = system('hostname');
disp(['My hostname is: ', result])
disp(['Attempting to perform nslookup on host: ', hostname])
disp('----------------------------------------------------------')
[status, result] = system(['nslookup ', hostname]);
disp(result);
disp('----------------------------------------------------------')

disp(['Attempting to perform ping on host: ', hostname])
disp('----------------------------------------------------------')
[status, result] = system(['ping ', hostname]);
disp(result);


portNumber = port;
disp('Attempting to connect back to the client using:')
fprintf('hostname: %s, portNumber: %d\n', hostname, portNumber);
disp('----------------------------------------------------------')
sa = java.net.InetSocketAddress(hostname, portNumber)
s = java.net.Socket()
try
    s.connect(sa, 120000)
catch err
    disp('Connection attempt failed');
    disp(err.getReport())
    try
        s.close();
    catch ME
        disp('Error closing the socket');
        disp(ME.getReport)
    end
    ok = false;
    disp('----------------------------------------------------------')
    
    return
end
disp('Connection attempt succeeded');
disp('----------------------------------------------------------')
try
    s.close();
catch ME
    disp('Error closing the socket');
    disp(ME.getReport)
end
ok = true;
return
end

% This function is run on the client. We establish a server socket and wait
% for the worker to connect back to us. The server socket has a timeout as
% it is a blocking operation.
function ok = establishServerSocket(port)
portNumber = port;
try
    s = java.net.ServerSocket(portNumber)
    s.setSoTimeout(120000) %3 minutes
    ns = s.accept()
catch err
    disp('Failed to establish a server socket on client')
    err.getReport()
    disp('----------------------------------------------------------')
    try
        ns.close()
    catch ME
        disp('Error closing the ns socket')
        disp(ME.getReport())
    end
    try
        s.close()
    catch ME
        disp('Error closing the s socket')
        disp(ME.getReport())
    end
    ok = false;
    return
end
try
    ns.getInetAddress().getHostAddress()
    ns.getPort()
catch err
    disp('Failed to obtain socket information on client')
    err.getReport()
    disp('----------------------------------------------------------')
    try
        ns.close()
    catch ME
        disp('Error closing the ns socket')
        disp(ME.getReport())
    end
    try
        s.close()
    catch ME
        disp('Error closing the s socket')
        disp(ME.getReport())
    end
    ok = false;
    return
end
try
    ns.close()
catch ME
    disp('Error closing the ns socket')
    disp(ME.getReport())
end
try
    s.close()
catch ME
    disp('Error closing the s socket')
    disp(ME.getReport())
end
ok = true;
disp('----------------------------------------------------------')
return
end

% iTestScheduler
% Determine whether the user specified configuration describes an actual
% scheduler.
function success = iTestParCluster(config)
try
    myCluster = parcluster(config)
    success = true;
catch err
    err.getReport('extended')
    success = false;
end
fprintf('\n')
end

% iTestDistributedJobs
% Test and verify that distributed jobs run properly
function success = iTestIndependentJobs(config)
timeOut = 600; %seconds
try
    myCluster = parcluster(config);
    job = createJob(myCluster);
    task = createTask(job, @rand, 1, {100});
    submit(job);
    ok = wait(job,'finished',timeOut);
    if ok && isempty(job.Tasks(1).ErrorIdentifier)
        success = true;
    else
        job.Tasks(1)
        success = false;
    end
catch err
    err.getReport('extended')
    success = false;
end
end

% iGatherLocalSystemInformation
% This function gathers client and worker operating system specific
% information that can be useful for debugging.
function iGatherLocalSystemInformation()
[hostnames, portrange] = gatherLocalHostNames();

% Obtain the ver of the system
disp('Obtaining ver')
disp('----------------------------------------------------------')
ver
fprintf('\n')
% Obtain information about MATLAB Path
disp('Obtaining MATLAB Path')
disp('----------------------------------------------------------')
path
fprintf('\n')

% Obtain information about startup.m and pathdef.m
disp('Obtaining information about startup.m')
disp('----------------------------------------------------------')
which -all startup.m
fprintf('\n')

disp('Obtaining information about pathdef.m')
disp('----------------------------------------------------------')
which -all pathdef.m
fprintf('\n')

% Determine if the system has gpuDevices
disp('Obtaining information about GPUs')
disp('This may take a few minutes')
disp('----------------------------------------------------------')
pause(1); %so everything updates on the screen
gpuSucceeded = false;
try
    disp('CUDA Driver Version')
    parallel.internal.gpu.CUDADriverVersion
    gpuSucceeded = true;
catch err
    disp('No GPU device drivers detected or error querying drivers')
    err.getReport('extended')
end

if gpuSucceeded
try
    disp('gpuDeviceCount')
    gpuDeviceCount
catch err
    disp('No GPU device drivers detected or error querying drivers/device')
    err.getReport('extended')
end

try
    disp('List all gpuDevices')
    for ii = 1:gpuDeviceCount
        fprintf('gpuDevice(%d)\n', ii);
        gpuDevice(ii)
    end
catch err
    disp('No GPU device drivers detected or error querying drivers/device')
    err.getReport('extended')
end
end

disp('WARNING: When using remote desktop or using workers managed by MJS')
disp('on Windows Vista/Windows 7, GPUs attached to display, and not running')
disp('TCC drivers will not be visible. (Sol # 1-FVMBP3)')
fprintf('\n')

% Determine additional system specific information
% Firewall settings
% All net processes and their associated ports
% IPCONFIG
% IFCONFIG
% Linux specific information
if(isunix)
    disp('This is a Linux or a Mac machine')
    % Obtain hostname
    disp('Obtaining hostname (hostname)')
    disp('----------------------------------------------------------')
    [~,b] = system('hostname');
    disp(b)
    % Obtaining system information
    disp('Obtaining system information (uname -a)')
    disp('----------------------------------------------------------')
    [~,b] = system('uname -a');
    disp(b);
    % Obtain environment
    disp('Obtaining environment')
    disp('----------------------------------------------------------')
    [~,b] = system('env');
    disp(b)
    % Obtain disk information
    disp('Obtaining disk information')
    disp('----------------------------------------------------------')
    [~,b] = system('df -h');
    disp(b)
    % Obtain system memory information
    disp('System memory information')
    disp('----------------------------------------------------------')
    [~,b] = system('free -m');
    disp(b);
    % Obtain information about CPU
    disp('Obtaining information about the CPUs')
    disp('----------------------------------------------------------')
    [~,b] = system('cat /proc/cpuinfo');
    disp(b);
    % Obtain information about processes running
    disp('Obtaining information about processes running on the system')
    disp('Output of top -n 1 -b')
    disp('----------------------------------------------------------')
    [~,b] = system('top -n 1 -b');
    disp(b)
    disp('Output of ps -ae')
    disp('----------------------------------------------------------')
    [~,b] = system('ps -ae'); %see which one displays full command info
    disp(b)
    % Obtain networking configuration
    disp('Obtaining networking information')
    disp('Obtaining contents of ifconfig')
    disp('----------------------------------------------------------')
    [~,b] = system('ifconfig');
    disp(b)
    
    disp('Obtaining contents of iptables')
    disp('----------------------------------------------------------')
    [~,b] = system('iptables -L -v');
    disp(b)
    % Obtain routing table
    disp('Obtaining routing table (netstat -r)')
    disp('----------------------------------------------------------')
    [~,b] = system('netstat -r');
    disp(b);
    % Obtain information about connection and associated programs
    disp('Obtaining information about network connections and associated programs (netstat -ab)')
    disp('----------------------------------------------------------')
    [~,b] = system('netstat -ab');
    disp(b);
    % Obtain protocol statistics
    disp('Obtaining protocol statistics (netstat -s)')
    disp('----------------------------------------------------------')
    [~,b] = system('netstat -s');
    disp(b);
    % Obtain system limitations
    disp('Obtaining system configuration and limits')
    disp('Output of sysctl -a')
    disp('----------------------------------------------------------')
    [~,b] = system('sysctl -a');
    disp(b);
    disp('Output of ulimit -a');
    disp('----------------------------------------------------------')
    [~,b] = system('ulimit -a');
    disp(b);
end
%% Windows specific information
if(ispc)
    disp('This is a windows machine')
    % Obtain hostname
    disp('Obtaining hostname (hostname)')
    disp('----------------------------------------------------------')
    [~,b] = system('hostname');
    disp(b)
    % Obtain environment
    disp('Obtaining environment')
    disp('----------------------------------------------------------')
    [~,b] = system('set');
    disp(b)
    % Obtain routing table
    disp('Obtaining routing table (netstat -r)')
    disp('----------------------------------------------------------')
    [~,b] = system('netstat -r');
    disp(b);
    % Obtain information about connection and associated programs
    disp('Obtaining information about network connections and associated programs (netstat -ab)')
    disp('----------------------------------------------------------')
    [~,b] = system('netstat -ab');
    disp(b);
    % Obtain protocol statistics
    disp('Obtaining protocol statistics (netstat -s)')
    disp('----------------------------------------------------------')
    [~,b] = system('netstat -s');
    disp(b);
    % Obtain network information
    disp('Obtaining network configuration (ipconfig)');
    disp('----------------------------------------------------------')
    [~,b] = system('ipconfig');
    disp(b)
    disp('Obtaining network configuration (netsh interface ip show config)');
    disp('----------------------------------------------------------')
    [~,b] = system('netsh interface ip show config');
    disp(b)
    disp('----------------------------------------------------------')
    
end
%% Mac specific information
if(ismac)
    disp('This is a Mac machine')
end

end
