function name = path2name(str, varargin)
%path2name(str...) take full path and return just distinctive name path

[a,b,c] = GetMonkeyName(str);
name = [a c];
