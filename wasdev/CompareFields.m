function [xa, xb] = CompareFields(a,b, varargin)
%  [uniquea, uniqueb] CompareFields(a,b, .)
%lists fields in a that are not in b, and vice versa

xa = {};
xb = {};
f = fields(a);
for j = 1:length(f)
    if ~isfield(b,f{j})
        xa = {xa{:} f{j}};
    end
end

f = fields(b);
for j = 1:length(f)
    if ~isfield(a,f{j})
        xb = {xb{:} f{j}};
    end
end
