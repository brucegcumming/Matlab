function [spike] = readNEV(nevfile)
%function [spike] = readNEV(nevfile)
%
% readNEV takes an NEV file as input and returns the event codes
% and times. 
%
% The columns of "spike" are channel, spike class (or digital port
% value for digital events), and time (seconds). The channel is 0
% for a digital event, and 1:255 for spike channels. 1:128 are the
% array channels, and 129 is the first analog input.
%
% This is a mex-optimized version of the matlab NEV_reader.m file
