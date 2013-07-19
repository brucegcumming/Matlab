function NEV = openNEV(NEV)

% Opens an .nev file for reading, returns all file information in a NEV
% structure. Works with File Spec 2.2.
% Use OUTPUT = openNEV(fname, 'read', 'report', 'electrodes', 'duration', 'mode').
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
%   'electrodes': User can specify which electrodes need to be read. The
%                 number of electrodes can be greater than or equal to 1
%                 and less than or equal to 128. The electrodes can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual electrodes (e.g. 3,6,7,90) or both.
%                 This field needs to be followed by the prefix 'e:'. See
%                 example for more details.
%                 DEFAULT: will read all existing channels.
%
%   'duration':   User can specify the beginning and end of the data
%                 segment to be read. If the start time is greater than the
%                 length of data the program will exit with an error
%                 message. If the end time is greater than the length of
%                 data the end packet will be selected for end of data. The
%                 user can specify the start and end values by comma 
%                 (e.g. [20,50]) or by a colon (e.g. [20:50]).
%                 This field needs to be followed by the prefix 't:'. See
%                 example for more details.
%                 DEFAULT: will read the entire file.
%
%   'mode':       The user can specify the mode of duration in [duration],
%                 such as 'sec', 'min', 'hour', or 'packet'. If 'sec' is
%                 specified the numbers in [duration] will correspond to
%                 the number of seconds. The same is true for 'min', 'hour'
%                 and 'packet'.
%                 DEFAULT: will be set to 'packet'.
%
%   OUTPUT:       Contains the NEV structure.
%
%   Example: 
%   
%   openNEV('report','read','c:\data\sample.nev', 'e:15:30', 't:3:10','min');
%
%   In the example above, the file c:\data\sample.nev will be used. A
%   report of the file contents will be shown. The data will be read from
%   electrodes 15 through 50 in the 3-10 minute time interval. If any of
%   the arguments above are omitted the default values will be used.
%
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 1.0.4 - March 09, 2009

%% Defining constants

[sortedStruct, sortedIndeces] = sortNEV(NEV.Data.SerialDigitalIO, 'Type', 'ExpParameter');
endIndeces = sortedIndeces - 1;
endIndeces(1) = [];

for i = 1:length(endIndeces)
    for j = 1:(endIndeces(i)-sortedIndeces(i)+1)
        organizedStruct{i,j} = NEV.Data.SerialDigitalIO(sortedIndeces(i)+j-1);
    end
end

NEV.Data.SerialDigitalIO = organizedStruct;

if ~nargout
    assignin('base', 'orgNEV', NEV);
    clear all;
end

%% Closing and clearing memory

end