function resps = runallsims(sims,varargin)
%resps = runallsims(sims,varargin)
% sims is a vector specifying which simulations to run. 
% resps is a cell array of results
% simulation 1 detect small V DG across all ori,dori 10k runs takes 
% simulation 2 detect small H DG across all ori,dsf
% simulation 3 look at dx/dg combos for V, H, across oris
% simulation 4 look at dx/dg combos for V, H, vertical ori
nrep = 10000;
saving = 1;
dg = 0.1;
j= 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nosave',2)
        saving = 0;
    elseif strncmpi(varargin{j},'nrep',2)
        j = j+1;
        nrep = varargin{j};
    end
    j = j+1;
end


for j = 1:length(sims)
    resps{j} = runsims(sims(j), nrep, saving);
end

function resps = runsims(sim, nrep, saving)
smalldg = 0.1;

tic;
if sim == 1
    resps = rundori('ori',[pi/16:pi/16:pi/2],'disps',0,'dgs',smalldg,'doris',[-pi/4:pi/16:pi/4],'nruns',nrep,'silent');
    resps.duration = toc;
    resps.runstart = now;
    save('oridori.mat','resps');
elseif sim == 2
    resps = rundori('ori',[pi/16:pi/16:pi/2],'disps',0,'dfs',smalldg,'dsf',[-1:0.1:1],'nruns',nrep,'silent');
    resps.duration = toc;
    resps.runstart = now;
    save('oridsf.mat','resps');
elseif sim == 3
    resps = rundori('ori',[pi/16:pi/16:pi/2],'disps',[-0.4:0.05:0.4],'nruns',nrep,'silent');
    resps.duration = toc;
    resps.runstart = now;
    save('ordxdor.mat','resps');
elseif sim == 4
    resps = rundori('ori',pi/2,'disps',[-0.4:0.05:0.4],'nruns',nrep,'track');
    resps.duration = toc;
    resps.runstart = now;
    if saving
    save('dxdor.mat','resps');
    end
end
