function str2 = LastConsecUpper(str1)
% returns the last occurence of 2 or more consecutive upper-case
% A-Z in str1
% NB has a bug so that if the last occurrence is followed by a
% single captial, it fails.
indx=find(str1>=65 & str1<=90);
j1 = max(find(diff(indx)>1))+1;
if isempty(j1)
    j1 = 1;
end
j2 = max(find(diff(indx)==1))+1;
str2 = str1([indx(j1):indx(j2)]);