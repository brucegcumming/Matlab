function NEV = openNEV(varargin)

% Opens an .nev file for reading, returns all file information in a NEV
% structure. Works with File Spec 2.2.
%
% Use OUTPUT = openNEV(fname, 'read', 'report', 'noparse').
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
%   'noparse':    The code will not parse the experimental parameters.
%                 See below for format.
%                 DEFAULT: will parse the parameters.
%
%   OUTPUT:       Contains the NEV structure.
%
%   Example: 
%   
%   openNEV('report','read','c:\data\sample.nev');
%

%   In the example above, the file c:\data\sample.nev will be used. A
%   report of the file contents will be shown. The digital data will be
%   parsed. The data needs to be in the proper format (refer below).
%
%   The experimental parameters need to be in the following format for the
%   code to properly parse them:
%
%   *Label:Parameter1=value1;Parameter2=value2;Parameter3=value3;#
%   *ExpParameter:Intensity=1.02;Duration=400;Trials=1;PageSegment=14;#
%
%   The above line is an "ExpParameter". The parameters are, "Intensity,
%   Duration, Trials, and PageSement." The values of those parameters are,
%   respectively, "1.02, 400, 1, and 14."
%
%   Another example:
%
%   *Marker:MarkerVal=10;#
%
%   The above line is a "Marker". The marker value is 10.
%
%   The label, parameter name, and values are flexible and can be anything.
%   The only required formatting is that the user needs to have a label and
%   a parameter name. The label is separated by a ':' and each parameter is
%   seperated from its values by a '='. The parameters are separated by
%   ';'. A '*' is required at the beginning of each segment and it needs
%   to end with a '#'.
%
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 3.0.0 - August 27, 2009

%% Defining structures
NEV = struct('MetaTags',[],'IOLabels',[], 'ElectrodesInfo', [], 'Data',[]);
NEV.MetaTags = struct('Subject', [], 'Experimenter', [], 'DateTime',[],...
    'SampleRes',[],'Comment',[],'FileTypeID',[],'Flags',[]);
NEV.Data = struct('SerialDigitalIO', [], 'Spikes', []);
NEV.Data.Spikes = struct('Timestamps',[],'Electrode',[],...
    'Unit',[],'Waveform',[]);
NEV.Data.SerialDigitalIO = struct('InputType',[],'TimeStamp',[],...
    'TimeStampSec',[],'Type',[],'Value',[]);

%% Validating input arguments
for i=1:length(varargin)
    if strcmpi(varargin{i}, 'report')
        Report = varargin{i};
    elseif strcmpi(varargin{i}, 'read')
        ReadData = varargin{i};
    elseif strcmpi(varargin{i}, 'noparse')
        ParseData = 'noparse';
    else
        temp = varargin{i};
        if length(temp)>3 && strcmpi(temp(end-3),'.')
            fname = varargin{i};
            if exist(fname, 'file') ~= 2
                display('The file does not exist.');
                clear all;
                return;
            end
        else
            if ~isnumeric(varargin{i})
                display(['Invalid argument ''' varargin{i} ''' .']);
            else display(['Invalid argument ''' num2str(varargin{i}) ''' .']);
            end
            clear all;
            return;
        end
    end
end

%% Defining and validating variables
if ~exist('fname', 'var')
    [fname, fpath] = uigetfile('D:\Data\*.nev');
else
    [fpath,fname, fext] = fileparts(fname);
    fname = [fname fext];
    fpath = [fpath '/'];
end
if ~exist('Report', 'var'); Report = 'noreport'; end
if ~exist('ReadData', 'var'); ReadData = 'noread'; end
if ~exist('ParseData', 'var'); ParseData = 'parse'; end
if fname==0; clear all; return; end;
tic;

%% Reading BasicHeader information from file
FID                       = fopen([fpath fname], 'r', 'ieee-le');
BasicHeader               = fread(FID, 336, '*uint8');
NEV.MetaTags.FileTypeID   = char(BasicHeader(1:8)');
NEV.MetaTags.Flags        = dec2bin(typecast(BasicHeader(11:12), 'uint16'),16);
fExtendedHeader           = double(typecast(BasicHeader(13:16), 'uint32'));
PacketBytes               = double(typecast(BasicHeader(17:20), 'uint32'));
TimeRes                   = double(typecast(BasicHeader(21:24), 'uint32'));
NEV.MetaTags.SampleRes    = typecast(BasicHeader(25:28), 'uint32');
t                         = typecast(BasicHeader(29:44), 'uint16');
NEV.MetaTags.Comment      = char(BasicHeader(77:332)');
ExtHeaderCount            = typecast(BasicHeader(333:336), 'uint32');
if strcmpi(NEV.MetaTags.FileTypeID, 'NEURALEV')
    disp('Filespec: 2.2.')
%    METATAGS = textread([fpath fname(1:end-8) '.sif'], '%s');
%    NEV.MetaTags.Subject      = METATAGS{3}(5:end-5);
%    NEV.MetaTags.Experimenter = [METATAGS{5}(8:end-8) ' ' METATAGS{6}(7:end-7)];
else
    disp('Filespec: 2.1.');
end

%% Parsing and validating FileSpec and DateTime variables
NEV.MetaTags.DateTime = [num2str(t(2)) '/'  num2str(t(4)) '/' num2str(t(1))...
    ' ' datestr(double(t(3)), 'dddd') ' ' num2str(t(5)) ':'  ...
    num2str(t(6)) ':'  num2str(t(7)) '.' num2str(t(8))] ;

%% Removing extra garbage characters from the Comment field.
NEV.MetaTags.Comment(find(NEV.MetaTags.Comment==0,1):end) = 0;

%% Recording after BasicHeader file position
fBasicHeader = ftell(FID); %#ok<NASGU>

%% Reading ExtendedHeader information
for ii=1:ExtHeaderCount
    ExtendedHeader = fread(FID, 32, '*uint8');
    PacketID = char(ExtendedHeader(1:8)');
    switch PacketID
        case 'ARRAYNME'
            NEV.ArrayInfo.ElectrodeName    = char(ExtendedHeader(9:end));
        case 'ECOMMENT'
            NEV.ArrayInfo.ArrayComment     = char(ExtendedHeader(9:end));
        case 'CCOMMENT'
            NEV.ArrayInfo.ArrayCommentCont = char(ExtendedHeader(9:end));
        case 'MAPFILE'
            NEV.ArrayInfo.MapFile          = char(ExtendedHeader(9:end));
        case 'NEUEVWAV'
            PacketID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.ElectrodeID     = PacketID;
            NEV.ElectrodesInfo{PacketID, 1}.ConnectorBank   = char(ExtendedHeader(11)+64);
            NEV.ElectrodesInfo{PacketID, 1}.ConnectorPin    = ExtendedHeader(12);
            NEV.ElectrodesInfo{PacketID, 1}.DigitalFactor   = typecast(ExtendedHeader(13:14),'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.EnergyThreshold = typecast(ExtendedHeader(15:16),'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.HighThreshold   = typecast(ExtendedHeader(17:18),'int16');
            NEV.ElectrodesInfo{PacketID, 1}.LowThreshold    = typecast(ExtendedHeader(19:20),'int16');
            NEV.ElectrodesInfo{PacketID, 1}.Units           = ExtendedHeader(21);
            NEV.ElectrodesInfo{PacketID, 1}.WaveformBytes   = ExtendedHeader(22);
        case 'NEUEVLBL'
            PacketID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.ElectrodeLabel = char(ExtendedHeader(11:26));
        case 'NEUEVFLT'
            PacketID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.HighFreqCorner = typecast(ExtendedHeader(11:14),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.HighFreqOrder  = typecast(ExtendedHeader(15:18),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.HighFilterType = typecast(ExtendedHeader(19:20),'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.LowFreqCorner  = typecast(ExtendedHeader(21:24),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.LowFreqOrder   = typecast(ExtendedHeader(25:28),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.LowFilterType  = typecast(ExtendedHeader(29:30),'uint16');
        case 'DIGLABEL'
            Label                                 = char(ExtendedHeader(9:24));
            Mode                                  = ExtendedHeader(25);
            NEV.IOLabels{Mode+1, 1} = Label;
        case 'NSASEXEV' %% Not implemented in the Cerebus firmware. 
                        %% Needs to be updated once implemented into the 
                        %% firmware by Blackrock Microsystems.
            NEV.NSAS.Freq          = typecast(ExtendedHeader(9:10),'uint16');
            NEV.NSAS.DigInputConf  = char(ExtendedHeader(11));
            NEV.NSAS.AnalCh1Conf   = char(ExtendedHeader(12));
            NEV.NSAS.AnalCh1Detect = typecast(ExtendedHeader(13:14),'uint16');
            NEV.NSAS.AnalCh2Conf   = char(ExtendedHeader(15));
            NEV.NSAS.AnalCh2Detect = typecast(ExtendedHeader(16:17),'uint16');
            NEV.NSAS.AnalCh3Conf   = char(ExtendedHeader(18));
            NEV.NSAS.AnalCh3Detect = typecast(ExtendedHeader(19:20),'uint16');
            NEV.NSAS.AnalCh4Conf   = char(ExtendedHeader(21));
            NEV.NSAS.AnalCh4Detect = typecast(ExtendedHeader(22:23),'uint16');
            NEV.NSAS.AnalCh5Conf   = char(ExtendedHeader(24));
            NEV.NSAS.AnalCh5Detect = typecast(ExtendedHeader(25:26),'uint16');
        otherwise
            display(['PacketID ' PacketID ' is invalid.']);
            clear all;
            return;
    end
end

%% Recording after ExtendedHeader file position and calculating Data Length
%  and number of data packets
fseek(FID, 0, 'eof');
fData = ftell(FID);
DataPacketCount = (fData - fExtendedHeader)/PacketBytes;
DataLen = PacketBytes - 8; %#ok<NASGU>

%% Reading data packets if 'read' is passed as an argument
fseek(FID, fExtendedHeader, 'bof');
Timestamps        = fread(FID, DataPacketCount, '*uint32', PacketBytes-4);
fseek(FID, fExtendedHeader+4, 'bof');
PacketIDs         = fread(FID, DataPacketCount, '*uint16', PacketBytes-2);
fseek(FID, fExtendedHeader+6, 'bof');
tempClassOrReason = fread(FID, DataPacketCount, '*uchar', PacketBytes-1);
fseek(FID, fExtendedHeader+8, 'bof');
tempDigiVals      = fread(FID, DataPacketCount, '*uchar', PacketBytes-1);
fseek(FID, fExtendedHeader+9, 'bof');
% tempThreshValues is not used because it has not been implemented in the
% firmware by Blackrock Microsystems.
tempThreshValues  = fread(FID, [5 DataPacketCount], '5*int16', PacketBytes-10); %#ok<NASGU>
fseek(FID, fExtendedHeader+8, 'bof');

%% Parse read digital data. Please refer to help to learn about the proper
% formatting if the data.
nonNeuralIndices  = find(PacketIDs == 0);
neuralIndices     = find(PacketIDs ~= 0);
nonNeuTimestamps  = Timestamps(nonNeuralIndices);
NeuTimestamps     = Timestamps(neuralIndices);
ElecNums          = PacketIDs(neuralIndices);
UnitNums          = tempClassOrReason(neuralIndices);
DigiValues        = char(tempDigiVals(nonNeuralIndices)');
if strcmp(ParseData, 'parse')
    AsteriskIndices   = find(DigiValues == '*');
    DataBegTimestamps = nonNeuTimestamps(AsteriskIndices);
    PoundIndices      = find(DigiValues == '#');
    ColonIndices      = find(DigiValues == ':');
    AfterAstAscii     = double(DigiValues(AsteriskIndices+1));
    NumAfterAstAscii  = find(AfterAstAscii >= 48 & AfterAstAscii <= 57);
    CharAfterAstAscii = find(AfterAstAscii >= 58);
    MarkerBegIndices  = AsteriskIndices(NumAfterAstAscii)+1;
    MarkerEndIndices  = PoundIndices(NumAfterAstAscii)-1;
    ParamsBegIndices  = AsteriskIndices(CharAfterAstAscii)+1;
    ParamsEndIndices  = PoundIndices(CharAfterAstAscii)-1;

    InsertionReason   = tempClassOrReason(find(PacketIDs == 0));
    Inputs = {'Digital'; 'AnCh1'; 'AnCh2'; 'AnCh3'; 'AnCh4'; 'AnCh5'; 'PerSamp'; 'Serial'};

    for i = 1:length(MarkerBegIndices)
        NEV.Data.SerialDigitalIO(NumAfterAstAscii(i)).Value(1,:) = DigiValues(MarkerBegIndices(i):MarkerEndIndices(i));
        NEV.Data.SerialDigitalIO(NumAfterAstAscii(i)).Type(1,:) = 'Marker';
    end
    for i = 1:length(ParamsBegIndices)-1
        Param        = DigiValues(ParamsBegIndices(i):ParamsEndIndices(i));
        NEV.Data.SerialDigitalIO(CharAfterAstAscii(i)).Value(1,:) = Param;
        NEV.Data.SerialDigitalIO(CharAfterAstAscii(i)).Type(1,:)  = DigiValues(ParamsBegIndices(i):ColonIndices(i)-1);
        SemiCIndices = find(NEV.Data.SerialDigitalIO(CharAfterAstAscii(i)).Value == ';');
        SemiCIndices = [find(NEV.Data.SerialDigitalIO(CharAfterAstAscii(i)).Value == ':') SemiCIndices];
        EqualIndices = find(NEV.Data.SerialDigitalIO(CharAfterAstAscii(i)).Value == '=');
        try
        for j = 1:length(EqualIndices)
            NEV.Data.SerialDigitalIO(CharAfterAstAscii(i)).(Param((SemiCIndices(j)+1):(EqualIndices(j)-1))) ...
                = str2num(Param((EqualIndices(j)+1):(SemiCIndices(j+1)-1)));
        end
        catch
            disp('There is an error in the formating of the digital data.');
            disp('Please refer to the help for more information on how to properly format the digital data for parsing.');
        end
    end
    % Populate the NEV structure with timestamps and inputtypes for the
    % digital data
%    c = num2cell(DataBegTimestamps); [NEV.Data.SerialDigitalIO(1:length(NEV.Data.SerialDigitalIO)).TimeStamp] = deal(c{1:end-1});
%    c = num2cell(DataBegTimestamps/NEV.MetaTags.SampleRes); [NEV.Data.SerialDigitalIO.TimeStampSec] = deal(c{1:end-1});
%    c = {Inputs{InsertionReason(AsteriskIndices)}}; [NEV.Data.SerialDigitalIO.InputType] = deal(c{1:end-1});
else
    NEV.Data.SerialDigitalIO.UnparsedData = DigiValues;
end

%% Populate the NEV structure with spike timestamps, electrode numbers
% and unit numbers
NEV.Data.Spikes.Timestamps = NeuTimestamps;
NEV.Data.Spikes.Electrode  = ElecNums;
NEV.Data.Spikes.Unit       = UnitNums;

%% Reads the waveforms if 'read' is passed to the function
if strcmp(ReadData, 'read')
    SpikeWaveform = fread(FID, [(PacketBytes-8)/2 DataPacketCount], ...
                    [num2str((PacketBytes-8)/2) '*int16'], 8)';
    SpikeWaveform(nonNeuralIndices, :) = [];
    NEV.Data.Spikes.Waveform           = SpikeWaveform;
end

%% Calculating the length of the data
fseek(FID, -PacketBytes, 'eof');
DataDuration = fread(FID, 1, 'uint32=>double');
%% Show a report if 'report' is passed as an argument
if strcmp(Report, 'report')
    
    disp( '*** FILE INFO **************************');
    disp(['File Path           = '  fpath]);
    disp(['File Name           = '  fname   ]);
    disp(['Data Duration (min) = ' num2str(round(DataDuration/NEV.MetaTags.SampleRes/60))]);
    disp(['Packet Counts       = ' num2str(DataPacketCount)]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');    
    disp(['File Type ID        = '         NEV.MetaTags.FileTypeID]);
    disp(['Sample Resolution   = ' num2str(NEV.MetaTags.SampleRes)]);
    disp(['Date and Time       = '         NEV.MetaTags.DateTime]);
    disp(['Comment             = '         NEV.MetaTags.Comment(1:64)   ]);
    disp(['                      '         NEV.MetaTags.Comment(65:128) ]);
    disp(['                      '         NEV.MetaTags.Comment(129:192)]);
    disp(['                      '         NEV.MetaTags.Comment(193:256)]);
end

if ~nargout
    assignin('base', 'NEV', NEV);
    fclose(FID);
    clear all;
else
    fclose(FID);
end

%% Display how fast the function was executed.
display(['The load time was ' num2str(toc, '%0.1f') ' seconds.']);
%% Closing and clearing memory

end