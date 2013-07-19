function data = NSX_read(NSX, chan, totchan, loc, time, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Greger Lab NSX binary file reader
% Revision 1.0
% NSX_read Read one or more concurrent channels from a neuroport recording
% NSX_read(NSX, chan, totchan, loc, time, size1, size2, load_time) 
% NSX = Structure of NSX file
% chan = first channel to read in
% totchan = total # of channels to read, reads them concurrently from the
% first channel
% loc = location to start read (default unit is miliseconds)
% time = amount of time(samples) to read in (default unit is miliseconds)
% size1 = (optional) Units for loc 
% size2 = (optional) Units for time
%  ---- Units - 's' for seconds, 'm' for minutes,'h' for hours, 'l' for miliseconds, and 'p'
%        for samples
% datatype = (optional) Load as double 'd', or 16-bit integers 'i', default
%             is double
% load_time = (optional) display load time
% Copyright Kyle Thomson - 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













if nargin > 5
    switch varargin{1}
        case 's'
            size_loc = 1; 
        case 'm'
            size_loc = 60;
        case 'h'
            size_loc = 3600;
        case 'p'
            size_loc = NSX.Period;
        otherwise
            size_loc = 1/1000;
    end;
else    
    size_loc = 1/1000;
end;

if nargin > 6
    switch varargin{2}
        case 's'
            size_time = 1;
        case 'm' 
            size_time = 60;      
        case 'h'
            size_time = 3600;
        case 'p'
            size_time = NSX.Period;
        otherwise
            size_time = 1/1000;
    end;
else
    size_time = 1/1000;
end;

if nargin > 7
    switch varargin{3}
        case 'i'
            datatype = 'short';
        case 'd'
            datatype = 'double';
        otherwise
            datatype = 'double';
    end;
else
    datatype = 'double';            
end;
        
loc = (1/NSX.Period)*size_loc*loc;
samples = (1/NSX.Period)*size_time*time;
tic;
STATUS = fseek(NSX.FID, NSX.Header+loc*NSX.Channel_Count*2+(chan-1)*2, 'bof');
if (STATUS == -1)
    ferror(NSX.FID)
    disp('Data Load Failed. No file associated.');       
    data = 'NaN';
    return
end;
ChanSkip = NSX.Channel_Count - totchan;
data = fread(NSX.FID,[totchan samples],[num2str(totchan) '*short=>' datatype],ChanSkip*2);
load_time = toc;
if nargin > 8
    disp(['Load time was ' num2str(load_time) ' seconds.']);
end;