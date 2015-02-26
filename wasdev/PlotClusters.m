function varargout = PlotClusters(varargin)
%PlotClusters(dir, ......)
%reads in Cluster files made by AllVPcs and plots up clusters for each
%expt/cell
%PlotClusters(dir, 'load') preloads all of the cluter files. Without this,
%files are loaded as needed.
%PlotClusters(dir, 'loadauto') also loads autoclusters
%PlotClusters(dir, 'loadautoonly') loads only autoclusters
%PlotClusters(dir, 'loadfullv') also forces loading of all the FullV files.
%This allows it to call AllVPcs quickly for any expt.  Not implemented for
%Utah arrays yet
%        BE SURE YOU HAVE ENOUGH MEMORY 


varargout = PC.PlotClusters(varargin{:});
