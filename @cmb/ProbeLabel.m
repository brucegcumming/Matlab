function str = ProbeLabel(Expt)
H = Expt.Header;
str = [];
if isfield(H,'cellnumber') && H.cellnumber > 0
str = sprintf('Cell%d (P%.1f)',H.cellnumber,H.probe);
elseif isfield(Expt.Header,'probe')
str = sprintf('P%.1f',H.probe);
end

