function size = StructSize(S)
%find out which elements of a structure take the space

f = fields(S);
X.(f{1}) = S.(f{1});
a = whos('X');
size(1) = a.bytes;
for j = 2:length(f)
    X = CopyFields(X,S,f(1:j));
    a = whos('X');
    size(j) = a.bytes;
end
    