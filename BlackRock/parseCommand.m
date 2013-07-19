% Parses out a string in the format specified below, into a structure.
%
% CommandType:ParamName1=ParamValue1;ParamName2=ParamValue2;
%
% The command string can have as many parameters in it as possible.
%
% Example:
%   ExpParameters:TaskType=Disc;StimPattern=0;SRad=10.0000;PRad=30.0000;Int
%   =0.2000;Xpos=0.0000;Ypos=5.0000;Dur=0.2000;
%
%   Kian Torab
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 1.2.1 - February 12, 2010

function [parsedCommand errorFlag] = parseCommand(inputString)

errorFlag = 0;
inputString(find(inputString == '#')) = [];
colonPos = find(inputString == ':', 1);
parsedCommand.Type = inputString(1:colonPos-1);
inputString = inputString(colonPos+1:end);
try
    splitString = regexp(inputString(1:end), ';', 'split');
    splitString(end) = [];
    splitParams = regexp(splitString', '=', 'split');
    splitParams = reshape([splitParams{:}], 2, length(splitString))';

    if isempty(splitString)
        parsedCommand = [];
        errorFlag = 1;
    else
        for idx = 1:length(splitString)
            [numVal, OK] = str2num(splitParams{idx, 2});
            if OK
                parsedCommand.(splitParams{idx,1}) = numVal;
            else
                parsedCommand.(splitParams{idx,1}) = splitParams{idx, 2};
            end
        end
    end
catch
    disp(['Cannot parse: ' inputString]);
    errorFlag = 1;
end