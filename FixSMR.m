function FixSMR(name)
%
% fixes errors in specifif matlab data files from Spike2

if strfind('ic242',name)
%ic242 was done before binoc sent pi * pursuitdir in expts that change pi
%so need to get sign of pi manually from actual change in FP.
    load('/data/icarus/242/ic242idx.mat');
    dy = Expt.Trials.dfy-Expt.Trials.fy;
    id = find(abs(Expt.Trials.pi) > 0.005);
    Expt.Trials.pi(id) = abs(Expt.Trials.pi(id)) .* sign(dy(id));
    save('/data/icarus/242/ic242idx.mat','Expt','Expts');
    
end