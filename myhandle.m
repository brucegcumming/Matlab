function true = myhandle(a)
%myhandle(a) Checks that a is a graphics handle AND is not zero
true = double(a) >0 & ishandle(a);
