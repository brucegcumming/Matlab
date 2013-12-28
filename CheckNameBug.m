function name = CheckNameBug(name)
%name = CheckNameBug(name)
%replaces leading '/ with '\' if ispc to undo the ridiculous bug in 2013b
%That is too hard for poor old matlb to sort out, especially consdidering
%the paltry sums paid to them by the US Government. 

if ispc && name(1) == '/'
    name(1) = '\';
end