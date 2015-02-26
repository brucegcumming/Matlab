function [type, details] = Type(name)
%Determines if name is online, Laminar probe files, Utah File, or standard
%Spike2

if iscellstr(name)
    for j = 1:length(name)
        [type{j}, details{j}] = expt.Type(name{j});
    end
    return;
end
type = '';
details.nfiles = [0 0 0];
if isdir(name)
    [a,b,c,d] = GetMonkeyName(name);
    d = dir([name '/*.mat']);
    nc = sum(CellToMat(regexp({d.name},'Expt.*ClusterTimes.mat')));
    ns = sum(CellToMat(regexp({d.name}, [a c '.[0-9]+.mat'])));
    no = sum(CellToMat(regexp({d.name}, [a c '.mat'])));
    if ns == 0 && no > 0 %old style
        type = 'spike2single'
    elseif ns > 2 && nc > ns/2 %offline by suffix
        type = 'laminar';
    elseif nc == 0 || nc < ns/4 %online, or single file
        type = 'online';
    end
    details.nfiles = [ns nc no];
end


