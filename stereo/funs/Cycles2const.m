function const = Cycles2const(cycles, filtersize)
%
%	freq = cycles2freq(cycles, filtersize)
%
%	cycles	:cycles of sinewave per filtersize
%
%	returns	:constant to be used in sin()/cos() functions

const = 2*pi*cycles/filtersize;
