function s = sprintcell(X, varargin)
%s = sprintcell(X, varargin) build a string from a cell array
% with multiple calls to sprintf;
s = [];
for j = 1:length(X)
    s = [s sprintf(varargin{:},X{j})];
end