function [nrow, ncol] = Nsubplots(nplots, varargin)
%[nrow, ncol] = Nsubplots(nplots, varargin)

if nplots > 36
    nrow = ceil(sqrt(nplots));
    ncol = ceil(nplots/nrow);
elseif nplots > 30
    nrow = 6;
    ncol = 6;
elseif nplots > 25
    nrow = 6;
    ncol = 5;
elseif nplots > 20
    nrow = 5;
    ncol = 5;
elseif nplots > 16
    nrow = 4;
    ncol = 5;
elseif nplots > 12
    nrow = 4;
    ncol = 4;
elseif nplots > 9
    nrow = 3;
    ncol = 4;
elseif nplots > 6
    nrow = 3;
    ncol = 3;
elseif nplots > 4
    nrow = 2;
    ncol = 3;
elseif nplots > 2
    nrow = 2;
    ncol = 2;
elseif nplots > 1
    nrow = 1;
    ncol = 2;
else
    nrow = 1;
    ncol = 1;
end
