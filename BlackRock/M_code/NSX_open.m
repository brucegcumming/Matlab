function [varargout] = NSX_open(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Greger Lab NSX binary file open handler
% Revision 1.0
% Open an NSX file for reading, returns all file information. 
% Currently works for File Spec 2.1 only. 
% Use NSX_open(fname, 1) to display file contents
% Copyright Kyle Thomson - 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%if no argument supplied, give GUI open
if nargin == 0
    [filename, pathname, filterindex] = uigetfile('*.ns*', 'Pick a Neural Recording');
    fname = [pathname filename];
else
    fname = varargin{1};
end

%Open as a MATLAB binary file
NSX.FID = fopen(fname,'r','ieee-le');

%Separate fname into path and filename, Windows and Unix compatible. 
backslash = find('\'== fname);
forwardslash = find('/'== fname);
if isempty(forwardslash)
    forwardslash = 0;
end;
if isempty(backslash)
    backslash = 0;
end;
Separator = max(backslash(end),forwardslash(end));
NSX.Filename = fname(Separator+1:end);
NSX.Path = fname(1:Separator);


%Read Binary file header
NSX.File_Type_ID = fread(NSX.FID, [1,8], 'char');
NSX.File_Spec = fread(NSX.FID, [1,16], 'char'); 
NSX.Period = 1/30000*fread(NSX.FID, 1, 'uint32');
NSX.Channel_Count = fread(NSX.FID, 1,  'uint32');

%Get the recording ID from each channel
for cnt = 1:NSX.Channel_Count
      NSX.Channel_ID(cnt) = fread(NSX.FID,1,'uint32');
end;

%Calculate header size
NSX.Header = 8+16+4+4+4*NSX.Channel_Count;

%Determine file length
fseek(NSX.FID,0,'eof');
NSX.Length = (ftell(NSX.FID) - NSX.Header)*NSX.Period / (NSX.Channel_Count*2);

%If no argument supplied, or 2, display file info. 
if (nargin > 1) || (nargin == 0)
    disp(['File Type ID = ' char(NSX.File_Type_ID)]);
    disp(['Label = ' char(NSX.File_Spec)]);
    disp(['Sampling Rate = ' num2str(1/NSX.Period)]);
    disp(['Channel Count = ' num2str(NSX.Channel_Count)]);
    disp(['Length = ' num2str(NSX.Length/60) ' minutes']);
end;

%If no output argument supplied, create a variable in the workspace equal
%to the file extension, otherwise, save to the supplied output argument. 
if nargout == 0
   varname = ['NS' fname(end)];
   assignin('base',varname,NSX);
else 
   varargout{1} = NSX;
end;

   

   
    