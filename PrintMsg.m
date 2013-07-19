function s = PrintMsg(logfid, varargin)
%s = PrintMsg(logfid,format, args)   calls sprintf(format,args)
%then prints the string to the console, and into the file handle logfid
%(if logfid > 0)
%returns printed message

s = sprintf(varargin{:});
fprintf('%s\n',s);
if logfid > 0
    try
    fprintf(logfid,'%s\r\n',s);
    catch
        fprintf('fid %d no longer valie\n',logfid);
    end
end