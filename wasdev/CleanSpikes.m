function ch = CleanSpikes(ch, varargin)
% CleanSpikes(ch, ...
% removes effects of saturation on spike channel in spike2 data files
% ch is a WaveMark channel as exported from SPike2 to matlab

bufl = 0;
profile = 0;
dvfile = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'buf',3)
        j = j+1;
        bufl = varargin{j};
    elseif strncmpi(varargin{j},'dvfile',3)
        j = j+1;
        dvfile = varargin{j};
    elseif strncmpi(varargin{j},'time',3)
        profile = 1;
    end
    j = j+1;
end
if profile
tic;
end
if dvfile & exist(dvfile,'file')
    fprintf('Loading %s...',dvfile);
    load(dvfile);
    v = dVdt;
elseif bufl > 0 & bufl < size(ch.values,1)
    bufpt = 1;
    lastpt = bufpt+bufl-1;
    while bufpt < size(ch.values,1)
    if lastpt > size(ch.values,1)
        lastpt = size(ch.values,1);
    end
    v(bufpt:lastpt,:,:) = diff(ch.values(bufpt:lastpt,:,:),1,2);
    bufpt = bufpt+bufl;
    lastpt = lastpt+bufl;
    end
else
    v = diff(ch.values,1,2);
end
[r,c] = find(abs(v) > 5); % change > 5v in one step = overflow
         % should add check that sign changes?
         ch.dVdt = v;
         clear v  %save some memory
if ndims(ch.dVdt) < 3
for spk = unique(r)'
    c = find(abs(ch.dVdt(spk,:)) > 5);
    for j = 1:2:length(c)
        if(length(c) == j) % last one never comes back
            c(j+1) = size(ch.values,2);
        end
        if ch.dVdt(spk,c(j)) < 0
            ch.values(spk,c(j)+1:c(j+1)) =  ch.values(spk,c(j)+1:c(j+1))+10;
        else
            ch.values(spk,c(j)+1:c(j+1)) =  ch.values(spk,c(j)+1:c(j+1))-10;
        end
    end
    ch.dVdt(spk,:) = diff(ch.values(spk,:));
end
end
if profile
toc
end
if dvfile & ~exist(dvfile)
    dVdt = ch.dVdt;
    fprintf('Saving dVdt to %s\n',dvfile);
    save(dvfile,'dVdt');
end