function Clusters = ReadTemplateClusters(a)
%ReadTemplateClustes(a)
%Reads a set of Cluster files build from templates by AllVPcs
%into a a struct array that so that set/get 'UserData' is not lsow

%with mode BAD, putting the structure into a figures UserDATA leads to very
%slow responses to get(f,'UserData') - c 20 secs!

GOOD = 1;
BAD = 2;
mode = GOOD;
F = {'templatexc' 'templatesrc' 'spts' 'probe' 'exptno' 'gmdprime' 'xy' 'savetime' 'Trigger' ...
    'bestspace' 'bestd' 'mahal' 'ctime' 'dropi'};

if mode == BAD
    F = {F{:} 'gmfit'};
end

for j = length(a):-1:1
    for k = length(a{j}):-1:1
        clear cluster;
        for c= length(a{j}{k}.cluster):-1:1
            cluster{c}.MeanSpike = a{j}{k}.cluster{c}.MeanSpike;
            for f = 1:length(F)
                if isfield(a{j}{k}.cluster{c},F{f})
                cluster{c}.(F{f}) = a{j}{k}.cluster{c}.(F{f});
                end
            end
        end
        Clusters{j}{k}.cluster = cluster;
    end
end
