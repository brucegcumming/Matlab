function NSx = saveNSx(NSx, varargin)

%%Opens an NSx file for reading, returns all file information in a NSx
% structure. Works with File Spec 2.1 and 2.2. 
% Use saveNSx(NSx, fname, 'report', 'channels', 'duration', 'mode').
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   NSx:          The data structure holding the channel information
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file. 
%                 DEFAULT: Will open Open File UI.
%
%   'report':     Will show a summary report if user passes this argument.
%                 DEFAULT: will not show report.
%
%   'channels': User can specify which channels need to be written. The
%                 number of channels can be greater than or equal to 1
%                 and less than or equal to 128. The channels can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual channels (e.g. 3,6,7,90) or both.
%                 This field needs to be followed by the prefix 'c:'. See
%                 example for more details.
%                 DEFAULT: will read all existing channels.
%
%   'duration':   User can specify the beginning and end of the data
%                 segment to be read. If the start time is greater than the
%                 length of data the program will exit with an error
%                 message. If the end time is greater than the length of
%                 data the end packet will be selected for end of data. The
%                 user can specify the start and end values by comma 
%                 (e.g. [20,50]) or by a colon (e.g. [20:50]). To use this
%                 argument the user must specify the [channels] or the
%                 interval will be used for [channels] automatically.
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
%   Example: 
%   
%   saveNSx(NSx, 'report','c:\data\sample.ns5', 'c:15:30', 't:3:10', 'min');
%   or equivalently
%   saveNSx(NSx, 'report','c:\data\sample.ns5', 'channels', 15:30, 'duration', 3:10, 'min');
%
%   In the example above, the file c:\data\sample.ns5 will be created. A
%   report of the file contents will be shown. The data of channels 15 
%   through 50 in the 3-10 minute time intervalwill be written to the file.
%   If any of the arguments above are omitted the default values will be used.
%
%   Original Author: Ehsan Azar
%
%   Contributors: 
%   Kian Torab, Blackrock Microsystems, kianabc@kianabc.com
%
%   Version 1.0.0.0
%

NSxver = '1.0.0.0';

%% Validating the input arguments. Exit with error message if error occurs.
for i=1:length(varargin)
    inputArgument = varargin{i};
    if strcmpi(inputArgument, 'report')
        Report = inputArgument;
    elseif strncmp(inputArgument, 't:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/'
        colonIndex = find(inputArgument(3:end) == ':');
        StartPacket = str2num(inputArgument(3:colonIndex+1));
        EndPacket = str2num(inputArgument(colonIndex+3:end));    
        if min(varargin{i})<1 || max(varargin{i})>128
            disp('The channel number cannot be less than 1 or greater than 128.');
            clear variables;
            if nargout; NSx = []; end
            return;
        end
    elseif strncmp(inputArgument, 'c:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/'
        Chan = str2num(inputArgument(3:end)); %#ok<ST2NM>        
    elseif strfind(' hour min sec sample ', [' ' inputArgument ' ']) ~= 0
        TimeScale = inputArgument;
    else
        temp = inputArgument;
        if length(temp)>3 && strcmpi(temp(end-3),'.')
            fname = inputArgument;
            if exist(fname, 'file') == 2
                disp('The file exists, it will be overwritten');
            end
        else
            disp(['Invalid argument ''' inputArgument ''' .']);
            clear variables;
            if nargout; NSx = []; end
            return;
        end
    end
end

tic;

%% Give all input arguments a default value. All input argumens are
%  optional.
if ~exist('Report', 'var');      Report = 'noreport'; end
if ~exist('StartPacket', 'var'); StartPacket = 0;     end
if ~exist('EndPacket', 'var'); EndPacket = length(NSx.Data);     end
if ~exist('TimeScale', 'var');   TimeScale = 'sample';       end
if ~exist('Chan', 'var');   Chan = 1:NSx.MetaTags.ChannelCount;       end

if (~isfield(NSx.MetaTags, 'Comment') || isempty(NSx.MetaTags.Comment))
    NSx.MetaTags.Comment = zeros(1, 256);
end

if (~isfield(NSx.MetaTags, 'DateTimeRaw') || isempty(NSx.MetaTags.DateTimeRaw))
    NSx.MetaTags.DateTimeRaw = zeros(8, 1);
end

if (~isfield(NSx.MetaTags, 'Timestamp') || isempty(NSx.MetaTags.Timestamp))
    NSx.MetaTags.Timestamp = 0;
end

% Defining constants
ExtHeaderLength            = 66;

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
    case 'sample'
        if (StartPacket > 0)
            StartPacket = StartPacket - 1;
        end
        EndPacket   = EndPacket - 1;
end

%% Validate StartPacket and EndPacket
if EndPacket >= length(NSx.Data)
    if StartPacket >= length(NSx.Data)
        disp('The starting packet is greater than the total data duration.');
        if nargout; NSx = []; end
        return;
    end
    disp('The time interval specified is longer than the data duration.');
    disp('Last data point will be used instead.');
    if strcmp(Report, 'report')
        disp('Press enter to continue...');
        pause;
    end
    EndPacket = length(NSx.Data) - 1;
else
DataLength = EndPacket - StartPacket + 1;

%% open file for write
if ~exist('fname', 'var')
    fname = [fullfile(NSx.MetaTags.FilePath, NSx.MetaTags.Filename) '.out.ns5'];
end

FID = fopen(fname, 'w', 'ieee-le');
if (FID <= 0)
    disp('Can not access file');
    return;
end

%% Writing Basic Header from NSx structure into file.
ChannelCount = length(Chan);

if (strcmp(NSx.MetaTags.FileSpec, '2.1') == 1)
    fwrite(FID, 'NEURALSG');
    fwrite(FID, NSx.MetaTags.SamplingLabel(1:16));
    fwrite(FID, NSx.MetaTags.TimeRes / NSx.MetaTags.SamplingFreq, 'uint32');
    fwrite(FID, ChannelCount, 'uint32');
    fwrite(FID, NSx.MetaTags.ChannelID(Chan), 'uint32');
else
    fHeader = 8 + 2 + 4 + 16 + 256 + 4 + 4 + 16 + 4 + ExtHeaderLength * ChannelCount;
    fwrite(FID, 'NEURALCD'); % 8
    fwrite(FID, [(NSx.MetaTags.FileSpec(1) - '0') (NSx.MetaTags.FileSpec(3) - '0')]); % 2
    fwrite(FID, fHeader, 'uint32'); % 4
    fwrite(FID, NSx.MetaTags.SamplingLabel(1:16)); %16
    if (length(NSx.MetaTags.Comment) < 256)
        NSx.MetaTags.Comment = [NSx.MetaTags.Comment zeros(256 - length(NSx.MetaTags.Comment), 1)];
    end
    strLen = find(NSx.MetaTags.Comment == 0, 1 , 'first');
    if (~isempty(strLen))
        NSx.MetaTags.Comment(strLen(1):end) = 0;
    end
    fwrite(FID, NSx.MetaTags.Comment(1:256)); % 256
    fwrite(FID, NSx.MetaTags.TimeRes / NSx.MetaTags.SamplingFreq, 'uint32'); % 4
    fwrite(FID, NSx.MetaTags.TimeRes, 'uint32'); % 4
    fwrite(FID, NSx.MetaTags.DateTimeRaw(1:8), 'uint16'); % 16 
    fwrite(FID, ChannelCount, 'uint32'); % 4
    % now write external header
    if (~isfield(NSx, 'ElectrodesInfo'))
        disp('Extended header information not found');
        fclose FID;
        return;
    end
    for ii = 1:ChannelCount
        fwrite(FID, 'CC'); % 2
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).ElectrodeID, 'uint16'); % 2
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).Label(1:16)); % 16
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).ConnectorBank - 'A' + 1); % 1
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).ConnectorPin); % 1
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).MinDigiValue, 'int16'); % 2
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).MaxDigiValue, 'int16'); % 2
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).MinAnalogValue, 'int16'); % 2
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).MaxAnalogValue, 'int16'); % 2
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).AnalogUnits(1:16)); % 16
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).HighFreqCorner, 'uint32'); % 4
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).HighFreqOrder,  'uint32'); % 4
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).HighFilterType, 'uint16'); % 2
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).LowFreqCorner, 'uint32'); % 4
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).LowFreqOrder,  'uint32'); % 4
        fwrite(FID, NSx.ElectrodesInfo(Chan(ii)).LowFilterType, 'uint16'); % 2
    end

    if fHeader ~= ftell(FID)
        disp('Header structure corrupted!');
        fclose(FID);
        return;
    end
    % data header
    fwrite(FID, 1); % 1
    fwrite(FID, NSx.MetaTags.Timestamp + StartPacket, 'uint32'); % 4
    fwrite(FID, DataLength, 'uint32'); % 4
end

%% write data points
fwrite(FID, NSx.Data(Chan, (StartPacket + 1):(EndPacket + 1)), 'int16');

%% Displaying a report of basic file information and the Basic Header.
if strcmp(Report, 'report')
    disp(['saveNSx ' NSxver]); clear NSxver;
    [path, fname, ~] = fileparts(fname);
    disp( '*** Input Data ***********************');
    disp(['Duration (seconds) = '  num2str(length(NSx.Data)/NSx.MetaTags.SamplingFreq)]);
    disp(['Total Data Points  = '  num2str(length(NSx.Data))]);
    disp(' ');
    disp( '*** FILE INFO **************************');
    disp(['File Path          = '  path]);
    disp(['File Name          = '  fname   ]);
    disp(['File Version       = '  NSx.MetaTags.FileSpec   ]);
    disp(['Duration (seconds) = '  num2str(DataLength/NSx.MetaTags.SamplingFreq)]);
    disp(['Data Points        = '  num2str(DataLength)]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');
    disp(['File Type ID       = '          NSx.MetaTags.FileTypeID      ]);
    disp(['Sampling Label     = '          NSx.MetaTags.SamplingLabel   ]);
    disp(['Sample Frequency   = '  num2str(NSx.MetaTags.SamplingFreq)         ]);
    disp(['Electrodes Write   = '  num2str(ChannelCount)   ]);
    
    disp(['The save time for ' fname ' file was ' num2str(toc, '%0.1f') ' seconds.']);
end

fclose(FID);

end