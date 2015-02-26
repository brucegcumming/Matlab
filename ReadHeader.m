function H = ReadHeader(name)

H = [];
if ~exist(name,'file')
    return;
end

a = matfile(name);
f = fields(a);
if sum(strcmp('SpikeHeader',f))
    H = a.SpikeHeader;
end