function x = CalcSubPlot(nr,nc, probe, Array, varargin)
%x = function CalcSubPlot(nr,nc, probe, Array)
%returns plot number to use for plotting probe given array, and
%number of rows/columns
j = 1;
while j <= length(varargin)
    j = j+1;
end


x = probe;
if ~isfield(Array,'X')
    return;
end

if nr == max(Array.X)
 x = Array.Y(probe) + nc * (nr-(Array.X(probe)));
end
