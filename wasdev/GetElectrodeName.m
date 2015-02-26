function estr = GetElectrodeName(Expt, varargin) 
%GetElectrodeName(Expt, varargin) find string describing electrode in Expt
estr = [];
if isfield(Expt,'Comments') && isfield(Expt.Comments,'Peninfo') && isfield(Expt.Comments.Peninfo,'electrode')
    estr = Expt.Comments.Peninfo.electrode;
end
