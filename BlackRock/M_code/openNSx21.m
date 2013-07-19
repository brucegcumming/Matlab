function NSx = openNSx(varargin)

%%Opens an NSx file for reading, returns all file information in a NSx
% structure. Works with File Spec 2.2. 
% Use OUTPUT = openNSx(fname, 'read', 'report', 'electrodes', 'duration', 'mode').
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file using an open 
%                 file user interface. 
%                 DEFAULT: Will open Open File UI.
%
%   'read':       Will read the data if user passes this argument.
%                 DEFAULT: will not read data.
%
%   'report':     Will show a summary report if user passes this argument.
%                 DEFAULT: will not show report.
%
%   [electrodes]: User can specify which electrodes need to be read. The
%                 number of electrodes can be greater than or equal to 1
%                 and less than or equal to 128. The electrodes can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual electrodes (e.g. 3,6,7,90) or both.
%                 This field needs to be followed by the prefix 'e:'. See
%                 example for more details.
%                 DEFAULT: will read all existing channels.
%
%   [duration]:   User can specify the beginning and end of the data
%                 segment to be read. If the start time is greater than the
%                 length of data the program will exit with an error
%                 message. If the end time is greater than the length of
%                 data the end packet will be selected for end of data. The
%                 user can specify the start and end values by comma 
%                 (e.g. [20,50]) or by a colon (e.g. [20:50]). To use this
%                 argument the user must specify the [electrodes] or the
%                 interval will be used for [electrodes] automatically.
%                 This field needs to be followed by the prefix 't:'. See
%                 example for more details.
%                 DEFAULT: will read the entire file.
%
%   'mode':       The user can specify the mode of duration in [duration],
%                 such as 'sec', 'min', 'hour', or 'sample'. If 'sec' is
%                 specified the numbers in [duration] will correspond to
%                 the number of seconds. The same is true for 'min', 'hour'
%                 and 'sample'.
%                 DEFAULT: will be set to 'sample'.
%
%   OUTPUT:       Contains the NSx structure.
%
%   Example: 
%   
%   openNSx('report','read','c:\data\sample.ns5', 'e:15:30', 't:3:10','min');
%
%   In the example above, the file c:\data\sample.ns5 will be used. A
%   report of the file contents will be shown. The data will be read from
%   electrodes 15 through 50 in the 3-10 minute time interval. If any of
%   the arguments above are omitted the default values will be used.
%
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 2.0.0 - June 1, 2009

%% Defining the NSx data structure and sub-branches.
NSx                = struct('MetaTags',[],'Data',[]);
NSx.MetaTags    = struct('FileTypeID',[],'ChannelCount',[],'SamplingFreq',[],'ChannelID',[]);

%% Validating the input arguments. Exit with error message if error occurs.
for i=1:length(varargin)
    inputArgument = varargin{i};
    if strcmpi(inputArgument, 'report')
        Report = inputArgument;
    elseif strcmpi(inputArgument, 'read')
        ReadData = inputArgument;
    elseif strncmp(varargin{i}, 't:', 2)
        inputArgument = str2num(inputArgument(3:end)); %#ok<ST2NM>
        StartPacket = inputArgument(1);
        EndPacket = inputArgument(end);      
        if min(varargin{i})<1 || max(varargin{i})>128
            display('The electrode number cannot be less than 1 or greater than 128.');
            clear all;
            return;
        end
    elseif strncmp(varargin{i}, 'e:', 2)
        Elec = str2num(inputArgument(3:end)); %#ok<ST2NM>        
    elseif strfind(' hour min sec sample ', [' ' inputArgument ' ']) ~= 0
        TimeScale = inputArgument;
    else
        temp = inputArgument;
        if length(temp)>3 && strcmpi(temp(end-3),'.')
            fname = inputArgument;
            if exist(fname, 'file') ~= 2
                display('The file does not exist.');
                clear all;
                return;
            end
        else
            display(['Invalid argument ''' inputArgument ''' .']);
            clear all;
            return;
        end
    end
end

%% Popup the Open File UI. Also, process the file name, path, and extension
%  for later use, and validate the entry.
if ~exist('fname', 'var')
    [fname, path] = uigetfile('D:\Data\*.ns*');
    if fname == 0
        clear all;
        return;
    end
    fext = fname(end-3:end);
else
    [path,fname, fext] = fileparts(fname);
    fname = [fname fext];
    path  = [path '\'];
end
if fname==0; return; end;

tic;

%% Give all input arguments a default value. All input argumens are
%  optional.
if ~exist('Report', 'var');      Report = 'noreport'; end
if ~exist('ReadData', 'var');    ReadData = 'noread'; end
if ~exist('StartPacket', 'var'); StartPacket = 0;     end
if ~exist('TimeScale', 'var');   TimeScale = 1;       end

%% Reading Basic Header from file into NSx structure.
FID                          = fopen([path fname], 'r', 'ieee-le'        );
NSx.MetaTags.FileTypeID   = fread(FID, [1,8]   , '*char'          );
Label                        = fread(FID, [1,16]  , '*char'          );
NSx.MetaTags.SamplingFreq = 30000/fread(FID, 1 , 'uint32=>double' );
ChannelCount                 = fread(FID, 1       , 'uint32=>double' );
NSx.MetaTags.ChannelID    = fread(FID, [ChannelCount 1], '*uint32');
NSx.MetaTags.ChannelCount = ChannelCount;
fHeader = ftell(FID);
fseek(FID, 0, 'eof');
fData = ftell(FID);
fseek(FID, fHeader, 'bof');

%% Validate the data file's File Spec. Exit with error message if invalid.
if ~strcmp(NSx.MetaTags.FileTypeID, 'NEURALSG')
    display('This version of openNSx can only read File Specs 2.1.');
    display(['The selected file label is ' NSx.MetaTags.FileTypeID '.']);
    fclose(FID);
    clear all;
    return; 
end;

%% Adjusts StartPacket and EndPacket based on what time setting (sec, min,
%  hour, or packets) the user has indicated in the input argument.
switch TimeScale
    case 'sec'
        StartPacket = StartPacket * NSx.MetaTags.SamplingFreq;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq;
    case 'min'
        StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 60;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 60;
    case 'hour'
        StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 3600;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 3600;
end

%% Validate StartPacket and EndPacket to make sure they do not exceed the
%  length of packets in the file. If EndPacket is over then the last packet
%  will be set for EndPacket. If StartPacket is over then will exist with an
%  error message.
NumofPackets = (fData-fHeader)/(2*ChannelCount);
if exist('EndPacket', 'var') && (EndPacket > NumofPackets)
    display('The time interval specified is longer than the data duration.');
    if StartPacket > NumofPackets
        clear all;
        return;
    end
    display('Last data point will be used instead.');
    display('Press enter to continue...');
    pause();
    EndPacket = NumofPackets;
elseif ~exist('EndPacket', 'var')
    EndPacket = NumofPackets;
end
DataLength = EndPacket - StartPacket;
clear TimeScale

%% Displaying a report of basic file information and the Basic Header.
if strcmp(Report, 'report')
    disp( '*** FILE INFO **************************');
    disp(['File Path          = '  path]);
    disp(['File Name          = '  fname   ]);
    disp(['Duration (seconds) = '  num2str(NumofPackets/NSx.MetaTags.SamplingFreq)]);
    disp(['Total Data Points  = '  num2str(NumofPackets)                   ]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');
    disp(['File Type ID       = '          NSx.MetaTags.FileTypeID      ]);
    disp(['Label              = '          NSx.MetaTags.Label           ]);
    disp(['Sample Resolution  = '  num2str(NSx.MetaTags.SamplingFreq)         ]);
    disp(['Electrodes Read    = '  num2str(NSx.MetaTags.ChannelCount)   ]);
end

%%
if ~exist('Elec', 'var'); Elec = 1:ChannelCount; end
fseek(FID, fHeader, 'bof');
if strcmp(ReadData, 'read')
    ReadElec = max(Elec)-min(Elec)+1;
    fseek(FID, StartPacket * 2 * ChannelCount, 'cof');
    fseek(FID, min(Elec)-1, 'cof');
    NSx.Data = fread(FID, [ReadElec DataLength-1], [num2str(ReadElec) '*int16'], ChannelCount-ReadElec);    
%     NSx.Data.DataPoints(end+1) = fread(FID, [ReadElec 1], '*int16');
end

%% If user does not specify an output argument it will automatically create
%  a structure.
if (nargout == 0),
    assignin('caller', ['NS' fext(4)], NSx);
end

display(['The load time was ' num2str(toc, '%0.1f') ' seconds.']);
fclose(FID);
clear all;

end