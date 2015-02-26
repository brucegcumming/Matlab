function probes = GetProbeList(DATA)
%return orderd list of which probes are being used, and
%order list to follow convention for TemplateScores,
% so if length >1, trigger channel is #2
% in returned list 1 is first non-trigger channel in list
%   2 is trigger channel
%   e is last non trigger channel in list
% 4-end any remaining channels

probes = DATA.chspk;
if isfield(DATA,'trueprobe')
    ispk = find(DATA.chspk == DATA.trueprobe);
else
    ispk = find(DATA.chspk == DATA.probe(1));
end

if isempty(ispk)
    ispk = 1;
end
if length(probes) > 1 && ispk == 1
    a = probes(2);
    probes(2) = probes(1);
    probes(1) = a;
end
if length(probes) > 2
    if ispk > 2
        a = probes(2);
        probes(2) = probes(ispk);
        probes(ispk) = a;
    end
else
    probes(3) = probes(end);
end
