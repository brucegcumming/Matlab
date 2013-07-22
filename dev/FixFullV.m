function FixFullV(path, varargin)

d = mydir([path '/Expt*.p1FullV.mat']);
for j = 1:length(d)
    load(d(j).name);
    if FullV.chspk == 96
        FullV.chspk = 1;
        save(d(j).name,'FullV');
    end
    clear FullV;
end
    
  