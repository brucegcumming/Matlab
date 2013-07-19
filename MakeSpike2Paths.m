function MakeSpike2Paths(name, netpref, varargin)
%Given pathname for a new Spike2 Data file, make sure necessary folders exist
%

[dname, fname] = fileparts(name);
[a, mname] = fileparts(dname);
monkey = GetMonkeyName(name);

odir = ['C:/Spike6/online/' monkey '/' mname];
fprintf('Creating %s,%s\n',odir,dname);
mkdir(odir);
mkdir(dname);
