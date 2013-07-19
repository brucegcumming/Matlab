function str2 = LastConsecUpper(str1)
% returns the last occurence of 2 or more consecutive upper-case
% A-Z in str1
indx=find(str1>=65 & str1<=90);
j1 = max(find(diff(indx)>1))+1;
j2 = max(find(diff(indx)==1))+1;
str2 = str1([indx(j1):indx(j2)]);