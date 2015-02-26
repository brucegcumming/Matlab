function result = DiskSummary(path)
%result = DiskSummary(path)
%searches path and summarized disk useage
%See also TreeSummary, which plots itself

[a,b,c,d] = TreeFind(path);
result.names = a;
result.sizes = b;
result.pathnames = d;
result.dates = c;

gb = (1024 * 1024 * 1024); 
totalbytes = sum(b);
fprintf('%.2f Gb\n',totalbytes./(1024 * 1024 * 1024));
types = {};
spkblk = [];
fullv = [];
for j = 1:length(a)
    if strfind(a{j},'.smr');
        types{j} = 'smr';
    elseif strfind(a{j},'FullV.mat');
        types{j} = 'FullV';
        id = strfind(a{j},'Expt');
        fullv(j) = sscanf(a{j}(id(1)+4:end),'%d');
    elseif strfind(a{j},'spkblk');
        types{j} = 'spkblk';
        id = regexp(a{j},'\.[0-9]*\.spkblk');
        spkblk(j) = sscanf(a{j}(id(1)+1:end),'%d');
    else
        types{j} = 'other';
    end
end
result.types = types;
alltypes = unique(types);
id = strmatch('spkblk',types);
spkblks = unique(spkblk(id));
id = strmatch('FullV',types);
fullvs = unique(fullv(id));

if length(spkblks)
    if length(fullvs)
        id = find(~ismember(spkblks,fullvs));
        if length(id)
            fprintf('Missing FullV: %s',sprintf(' %d',spkblks(id)));
        end
    else
        fprintf('No FullV Files\n');
    end
end

for j = 1:length(alltypes);
    id = strmatch(alltypes{j},types);
    sz = sum(b(id))./gb;
    if strcmp(alltypes{j},'FullV')
        fprintf('%s: %d files (%d-%d), %.3fGb\n',alltypes{j},length(id),min(fullvs),max(fullvs),sz);
    elseif strcmp(alltypes{j},'spkblk')
        fprintf('%s: %d files, %d expts (%d-%d), %.3fGb\n',alltypes{j},length(id),length(spkblks),min(spkblks),max(spkblks),sz);
    else
    fprintf('%s: %d files, %.3fGb\n',alltypes{j},length(id),sz);
    end
end
