 function val = GetTextVal(varargin)
% GetTextVal(args) used findobj(args) to locat a text GUI field,
% e.g. GetTextVal(F, 'Tag','smooth')  finds the child of F with Tag
% 'smooth', and returns the text converted to a number
% read the string, and converts to double
% returns NaN if no object is found


     val = NaN;
     it = findobj(varargin{:});
     for j = 1:length(it)
         str = get(it(j),'string');
         val(j) = str2num(str);
     end