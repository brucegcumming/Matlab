function [allframes] = RlsRCa(Expt, seedlist, idlist, varargin)
%
% allframes = RlsRCa(Expt, seedlist, idlist, varargin)
% Finds members of seedlist/idlist that match spikes in Expt, for
% various delays ([40 50 60]ms by default).
%
% The return value is a cell array, each cell containing the idlist for a
% given delay
% RlsRCa(Expt, seedlist, idlist, 'delays',X)
% uses the elements of vector X as the set of delays. 
%
% RlsRCa(Expt, seedlist, idlist, 'silence',X)
% returns all frame associated with a period of X tics with no spikes
% N.B. if the silence lasts longer than X, each frame associated with the
% additional silent period is returned. ...'silence',X,'exact')
% returns only the first frame

checkimages = 0;
firststim = 1;
laststim = length(Expt.Trials);
sumall = 0; %% build mean image, regardless of spikes
delays = [400 500 600];
prefix = '';
dotw = 10;
scale = 2;
silence = 0;
seedstep = 2; % for uncorrelated, need 1 for correlated
uselastseed = 0;
seedoffset = 0;
exact = 0;

j = 1;
while j < nargin-2
    if strncmpi(varargin{j},'check',4)
        checkimages = 1;
        if isnumeric(varargin{j+1})
            j = j+1;
            firststim = varargin{j};
            laststim = firststim+1;
        end
    elseif strncmpi(varargin{j},'all',3)
        sumall = 1;
    elseif strncmpi(varargin{j},'delays',4)
        j = j+1;
        delays = varargin{j};
    elseif strncmpi(varargin{j},'dotw',4)
        j = j+1;
        dotw = varargin{j};
    elseif strncmpi(varargin{j},'exact',4)
        exact = 1;
    elseif strncmpi(varargin{j},'step',4)
        j = j+1;
        seedstep = varargin{j};
    elseif strncmpi(varargin{j},'offset',4)
        j = j+1;
        seedoffset = varargin{j};
    elseif strncmpi(varargin{j},'prefix',4)
        j = j+1;
        prefix = varargin{j};
    elseif strncmpi(varargin{j},'silence',4)
        j = j+1;
        silence = varargin{j};
    elseif strncmpi(varargin{j},'forcelast',6)
        uselastseed = 2;
    elseif strncmpi(varargin{j},'uselast',5)
        uselastseed = 1;
    end
    j = j+1;
end

%imload('/bgc/bgc/anal/rc/duf1172im.mat')

frameperiod = 10000 / Expt.Stimvals.fz;
sframes = silence/frameperiod;
totalspikes = zeros(size(delays));
nullim = totalspikes;
nspikes = sum([Expt.Trials.Count]);
nallim = 0;
nframes = 0;
frame = 0;
dw = fix(dotw/scale);
kernel = ones(dw,dw)./(dw * dw);
%since dot is in fact finite size, adjacent dot locations all illuminate a
%pixel. So if dw is real dot size, don't need to divide kernel
kernel = ones(dw,dw);
if length(delays) > 6
    nrow = 3;
    ncol = 3;
elseif length(delays) > 4
    nrow = 2;
    ncol = 3;
else
    nrow = 2;
    ncol = 2;
end
usedseeds = [];
unused = [];

for j = 1:length(delays)
    allframes{j} = [];
%    allframes{j}(nspikes) = NaN; %% Allocate memory
end

for stim = firststim:laststim
    if ~isempty(Expt.Trials(stim).id)
        sid = Expt.Trials(stim).id;
%
% setting uselastseed to 2 forces the use of ls regardless, to deal with
% duf1191 where ls = se, but ls is correct;
        if uselastseed == 2 | (uselastseed & Expt.Trials(stim).ls > Expt.Trials(stim).se + Expt.Trials(stim).Nf);
            seed = Expt.Trials(stim).ls - Expt.Trials(stim).Nf .* seedstep;
        else
            seed = Expt.Trials(stim).se + seedoffset;
        end
        startseed = seed;
        if silence
            isi = diff(Expt.Trials(stim).Spikes);
            isid = find(isi > silence);
            spks = [];
            for k = 1:length(delays)
                seeds = [];
                for j = 1:length(isid)
                    sframes = floor((isi(isid(j)) - silence)./frameperiod);
                    frames = (Expt.Trials(stim).Spikes(isid(j)) + silence - delays(k))/frameperiod;
                    if ~exact
                        frames = frames + [0:sframes];
                    end
                    seeds = [seeds startseed + floor(frames) * seedstep];
                    %                spks = [spks (Expt.Trials(stim).Spikes(isid(j))+[silence:frameperiod:isi(isid(j))])];
                end
                ids  = find(ismember(seedlist,seeds));
                allframes{k} = [allframes{k} ids'];
            end
        else
            spks = Expt.Trials(stim).Spikes;
            seeds = [];
            idid = find(idlist == sid);
            idid = reshape(idid,length(idid),1);
            frames = spks./frameperiod;
            tseed = seedlist(idid)+seedoffset;
            for j = 1:length(delays)
                frames = (spks - delays(j))./frameperiod;
                seeds = startseed + floor(frames(find(frames > 0))) * seedstep;
                [a, id] = ismember(seeds, tseed);
                allframes{j} = [allframes{j} idid(id(find(id)))'];
            end
        end
    end
end


function [pixels, seed] = GetPixels(imname, imdim)

load(imname);
seed = RdsDots.seed;
pixels(:,1) = RdsDots.left.c;
pixels(:,2) = RdsDots.right.c;