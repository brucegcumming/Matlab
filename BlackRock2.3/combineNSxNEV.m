%
%   [NSx, NEV] = combineNSxNEV(filename1, filename2)
%
%   This function loads two NSx and NEV files and it will combine them
%   together into one file. The resulting file will be saved as new NSx
%   and NEV files.
%
%
%   filename1:  The name of the first NSx file. This input is optional. In
%               its absense, a dialog will open and will prompt the user to
%               select an NSx file.
%               (OPTIONAL)
%
%   filename2:  The name of the second NSx file. This input is also
%               optional. In its absense, a dialog will open and will
%               prompt the user to select an NSx file.
%               (OPTIONAL)
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%
%   Version 1.0.0.0


function [NSx1, NEV1] = combineNSxNEV(filename1, filename2)

% Openning NSx files
if exist('filename1', 'var') && exist('filename2', 'var')
    NSx1 = openNSx('read', filename1);
    NSx2 = openNSx('read', filename2);
    if NSx1.MetaTags.SamplingFreq ~= NSx2.MetaTags.SamplingFreq;
        disp('The sampling frequencies are not the same.');
        return;
    end
else
    NSx1 = openNSx('read');
    NSx2 = openNSx('read');
end

% Determining length of the first NSx file
conversionFactor = 30000/NSx1.MetaTags.SamplingFreq;
NSx1DataLength = NSx1.MetaTags.DataPoints * conversionFactor;

% Combining NSx files
NSx1.Data = [NSx1.Data, NSx2.Data];

% Opening NEV files
fileNameNEV1 = [NSx1.MetaTags.FilePath '/' NSx1.MetaTags.Filename(1:end-3), 'nev'];
fileNameNEV2 = [NSx1.MetaTags.FilePath '/' NSx2.MetaTags.Filename(1:end-3), 'nev'];
clear NSx2;

if (exist(fileNameNEV1, 'file') == 2) && (exist(fileNameNEV2, 'file') ==2)
    disp('Openning corresponding NEV files...');
    NEV1 = openNEV('read', fileNameNEV1);
    NEV2 = openNEV('read', fileNameNEV2);
else
    NEV1 = openNEV('read');
    NEV2 = openNEV('read');
end

% Adjusting the timestamp on the second NEV file
NEV2.Data.Comments.TimeStamp = NEV2.Data.Comments.TimeStamp + NSx1DataLength;
NEV2.Data.Comments.TimeStampSec = NEV2.Data.Comments.TimeStampSec + double(NSx1DataLength)/30;
NEV2.Data.SerialDigitalIO.TimeStamp = NEV2.Data.SerialDigitalIO.TimeStamp + NSx1DataLength;
NEV2.Data.SerialDigitalIO.TimeStampSec = NEV2.Data.SerialDigitalIO.TimeStampSec + double(NSx1DataLength)/30;
NEV2.Data.Spikes.TimeStamp = NEV2.Data.Spikes.TimeStamp + NSx1DataLength;
NEV2.Data.VideoSync.TimeStamp = NEV2.Data.VideoSync.TimeStamp + NSx1DataLength;
NEV2.Data.Tracking.TimeStamp = NEV2.Data.Tracking.TimeStamp + NSx1DataLength;
NEV2.Data.PatientTrigger.TimeStamp = NEV2.Data.PatientTrigger.TimeStamp + NSx1DataLength;
NEV2.Data.Reconfig.TimeStamp = NEV2.Data.Reconfig.TimeStamp + NSx1DataLength;

% Combining the two NEV files
NEV1.Data.Spikes.Electrode      = [NEV1.Data.Spikes.Electrode, NEV2.Data.Spikes.Electrode];
NEV1.Data.Spikes.TimeStamp      = [NEV1.Data.Spikes.TimeStamp, NEV2.Data.Spikes.TimeStamp];
NEV1.Data.Spikes.Unit           = [NEV1.Data.Spikes.Unit, NEV2.Data.Spikes.Unit];
NEV1.Data.Spikes.Waveform       = [NEV1.Data.Spikes.Waveform, NEV2.Data.Spikes.Waveform];
NEV1.Data.Comments.TimeStamp    = [NEV1.Data.Comments.TimeStamp, NEV2.Data.Comments.TimeStamp];
NEV1.Data.Comments.TimeStampSec = [NEV1.Data.Comments.TimeStampSec, NEV2.Data.Comments.TimeStampSec];
NEV1.Data.Comments.CharSet      = [NEV1.Data.Comments.CharSet, NEV2.Data.Comments.CharSet];
NEV1.Data.Comments.Color        = [NEV1.Data.Comments.Color, NEV2.Data.Comments.Color];
NEV1.Data.Comments.Text         = [NEV1.Data.Comments.Text; NEV2.Data.Comments.Text];
NEV1.Data.Tracking.TimeStamp    = [NEV1.Data.Tracking.TimeStamp, NEV2.Data.Tracking.TimeStamp];
NEV1.Data.Tracking.PointCount   = [NEV1.Data.Tracking.PointCount, NEV2.Data.Tracking.PointCount];
NEV1.Data.Tracking.X            = [NEV1.Data.Tracking.X, NEV2.Data.Tracking.X];
NEV1.Data.Tracking.Y            = [NEV1.Data.Tracking.Y, NEV2.Data.Tracking.Y];
NEV1.Data.Tracking.TimeStampSec = [NEV1.Data.Tracking.TimeStampSec, NEV2.Data.Tracking.TimeStampSec];
