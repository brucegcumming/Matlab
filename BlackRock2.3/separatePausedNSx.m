function separatePausedNSx(varargin)

% separatePausedNsx    Saves paused data segments from a single NSx file 
%                      as individual NSx files
%
%    separatePausedNSx(FILENAME) where FILENAME has N paused segments
%    will create N individual files with the same extension named
%    FILENAME-1, FILENAME-2, ..., FILENAME-N.
%
%    separatePausedFiles without any input arguments opens a UIgetfile to
%    select the NSx file to separate
%
%    Brett Dowden
%    bdowden@blackrockmicro.com
%    Blackrock Microsystems
%    Version 1.0.1.0

% Since openNSx checks input parameters for file validity and for no input
% argument case, just pass varargin to openNSx
NSx = openNSx('read',varargin{:});

% Check to see if there were pauses in the file.  If so, save the sections
% as individual files, if not just report that to the user and end function
if length(NSx.MetaTags.Timestamp) == 1
    disp('No pauses in data file.  No action taken');
else
    for i = 1 : length(NSx.MetaTags.Timestamp)
        % create a copy of NSx, except with only one data segment and the
        % corresponding Timestamp, DataPoints, and NumofPacket fields
        NSx_out.MetaTags              = NSx.MetaTags;
        NSx_out.MetaTags.Timestamp    = NSx.MetaTags.Timestamp(i);
        NSx_out.MetaTags.DataPoints   = NSx.MetaTags.DataPoints(i);
        NSx_out.MetaTags.NumofPackets = NSx.MetaTags.DataPoints(i);
        
        NSx_out.ElectrodesInfo        = NSx.ElectrodesInfo;
        
        NSx_out.Data                  = NSx.Data{i};
        
        % Create enumerated file name and save this segment as its own file
        NSx_out.MetaTags.Filename     = [NSx.MetaTags.Filename(1:end-4) '-p' sprintf('%03d', i) NSx.MetaTags.FileExt];
        saveNSx(NSx_out,NSx_out.MetaTags.Filename);
    end
end