function varargout = argtest(varargin)

if nargout == 1
varargout{1} = testarg(varargin{:});
elseif nargout == 2
[varargout{1}, varargout{2}] = testarg(varargin{:});
end


function varargout = testarg(varargin)

nargout
varargout{1} = {'test'};
varargout{2} = {'test2'};
nargout


