function loaddetails
% Illustrating another 2013b bug. this time a variable 'details'
% gets masked by some matlab script. But only if the variable
% details arose from loading a file. Put at stop at last line 
% copy and paste to command window works fine.  Then step over
% - error!!!!!!!!!!!!
%
%This code is fine in 2012b

%details.x = 1; %uncomment this line and it works!!
MakeDetails('Dvariable');
load('Dvariable','details'); %file containing the same details variable
who
which details
Kernels{1} = details;



function MakeDetails(name)

details.x = 1;
save(name,'details');