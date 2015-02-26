function runtest(name, varargin)

DATA.test = 1;
DATA.name = name;
DATA.clustersubdir = [];
if length(varargin) > 0 && isnumeric(varargin{1}) && varargin{1} > 0
    OTest.runtest(name, varargin{1}-1);
end
OTest.ReadSpikeFiles(DATA, name );
