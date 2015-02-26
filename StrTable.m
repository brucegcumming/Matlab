function T = StrTable(s, varargin)
%print a cell array of strings as a table
spacing = 1;

for j = 1:size(s,1)
    for k = 1:size(s,2)
        sizes(j,k) = length(s{j,k});
    end
end

flen = max(sizes) + spacing;

for j = 1:size(s,1)
    for k = 1:size(s,2)
        pad = [];
        n = flen(k) - sizes(j,k);
        pad = repmat(' ', n, 1);
        fprintf('%s%s',s{j,k},pad);
    end
    fprintf('\n');
end
