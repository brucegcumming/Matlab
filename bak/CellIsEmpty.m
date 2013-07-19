function true = CellIsEmpty(S, field, n)
%true CellIsEmpty(S, field, n) Checks struct S to see if cell array
%field has this element alreay. Test 1 field exits, 2 field lengtt, 3
%filend{n} is empty.

if ~isfield(S,field) || length(S.(field)) < n || isempty(S.(field){n})
    true = 1; %Cell is empty
else
    true = 0;
end