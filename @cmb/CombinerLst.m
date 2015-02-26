function cfile = CombinerLst(DATA)
if isfield(DATA,'datafilename') %%not online data
cfile = strrep(DATA.datafilename,'.mat','.combine.mat');
else
cfile = [DATA.datafilename '/combine.mat'];
end

