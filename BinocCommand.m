function BinocCommand(str, varargin)
% BinocCommand(str) Sends a single commnad to binoclean
%useful to co-ordinate gamma calibration

DATA.network = 1;
DATA.verbose = 1;
DATA.ip = 'http://localhost:1110/';
DATA.outpipe = '/tmp/binocinputpipe';
DATA.inpipe = '/tmp/binocoutputpipe';

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'network',5)
        DATA.network = 1;
    elseif strncmpi(varargin{j},'pipes',4)
        DATA.network = 0;
    end
    j = j+1;
end

if DATA.network == 0
    DATA.outid = fopen(DATA.outpipe,'w');
end
if iscellstr(str)
    for j = 1:length(str)
        outprintf(DATA,str{j});
    end
else
    outprintf(DATA,str);
end
 



function outprintf(DATA,varargin)
if DATA.network == 0
    pipeprintf(DATA,varargin{:});
    return;
end
 %send to binoc.  ? reomve comments  like in expt read? 
     str = sprintf(varargin{:});
     strs = split(str,'\n');
     for j = 1:length(strs)
         if ~isempty(strs{j})
             str = [DATA.ip strs{j}];
             if DATA.verbose
                 fprintf('%s\n',str);
             end
             ts = now;
             [bstr, status] = urlread(str,'Timeout',2);
             if ~isempty(bstr)
                 fprintf('Binoc replied with %s\n',bstr);
             elseif DATA.verbose
                 fprintf('Binoc returned in %.3f\n',mytoc(ts));
             end
         end
     end
     
function pipeprintf(DATA,varargin)
 %send to binoc.  ? reomve comments  like in expt read? 
     str = sprintf(varargin{:});
     strs = split(str,'\n');
     for j = 1:length(strs)
         if ~isempty(strs{j})
             fprintf(DATA.outid,'%s\n',strs{j});
         end
     end
