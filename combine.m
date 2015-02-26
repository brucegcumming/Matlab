function varargout = combine(varargin)
%cmb.combine(file)
%cuts clusters and combines expts for Spike2 generated matlab data files
%
%

% can't get rid of cluster 2 (M019)
% show cluster colors in time display
% build voltage average for synchronous spikes;
% add checkbox for plotflip
% deleting clusters in AllClusters box
%
% Show spikes for non-stim trials only
%right button on a rotate ellipase sometimes cancels angle
% setting cluster 2 in a new space erates cluster 1
%scaling of PCA often mucked up
%artifact setting second regios sometimes undoes first, lemM102 p3 ex19
%autoscale density plot to peak of cluster? 

if nargout == 0
    cmb.combine(varargin{:});
elseif nargout  == 1
    varargout{1} = cmb.combine(varargin{:});
elseif nargout  == 2
    [varargout{1}, varargout{2}] = cmb.combine(varargin{:});
end

