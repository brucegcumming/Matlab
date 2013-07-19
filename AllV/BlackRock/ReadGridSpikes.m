function DATA = ReadGridSpikes(DATA, name, varargin)

%ReadGridSpikes(DATA, name, varargin) 
%Reads a BlckRock NEV file into the DATA.AllSpikes structure
verbose = 0;
plotsummary = 0;
toff = 2577400;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plotspike',8)
        plotsummary =1;
    elseif strncmpi(varargin{j},'toff',4)
        j = j+1;
        toff = varargin{j};
    end
    j = j+1;
end
nsfile = strrep(name,'.nev','.ns5');
if exist(nsfile,'file')
    nev = ns2nev(nsfile);
else
    nev = openNEV('read',name);
end
DATA.AllSpikes = {};
np = floor(max(nev.Data.Spikes.Electrode));
np = min([np 96]); %kludge for now
trate = nev.MetaTags.TimeRes;
for j = 1:np
    id = find(nev.Data.Spikes.Electrode == j);
    if length(id)
        DATA.AllSpikes{j}.values = nev.Data.Spikes.Waveform(id,:)./512;
        DATA.AllSpikes{j}.codes = zeros(length(id),4);
        DATA.AllSpikes{j}.codes(:,1) = nev.Data.Spikes.Unit(id);
        DATA.AllSpikes{j}.times = nev.Data.Spikes.TimeStamp(id).* (10000/trate) + toff;
        DATA.AllSpikes{j}.firstspki = 1;
        DATA.Spikes.minv(j) = min(DATA.AllSpikes{j}.values(:));
    end
    DATA.nspk(j) = length(id);
end
if verbose
end
if plotsummary
    p = 0;
    GetFigure('GridSpikes');
    for j = 1:np
        p = p+1;
        subplot(10,10,p);
        if isempty(DATA.AllSpikes{j})
            plot(0,0);
        else
        v = min(DATA.AllSpikes{j}.values,[],2);
        id = find(v < prctile(v,10));
        if length(id)
        plot(DATA.AllSpikes{j}.values(id,:)','color',[0.5 0.5 0.5]);
        title(num2str(length(id)));
        else
            plot(0,0,'.');
        end
        end
        set(gca,'xticklabel',[]);
    end
end
