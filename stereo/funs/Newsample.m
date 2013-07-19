function im=Newsample(n)
%
% im=newsample(n)
%
% Makes a random dot image
%   basically returns 1d array
%   of contrast -1:+1
dim = [1 n];
im=(rand(dim));
im=round(im);
im = im * 2 -1;

