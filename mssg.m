function s = mssg(fid, varargin)

s = sprintf(varargin{:});
fprintf([s '\n']);
if fid > 0
    fprintf(fid,s);
end
