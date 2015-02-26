function x = myround(x,sd)

% allows you to specify the number of decimal places to round to

% Adrian Bondy, 2012

x=x.*10^sd;
x=round(x);
x=x./10^sd;
