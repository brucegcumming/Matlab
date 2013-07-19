function Expt = ReadOnline(name, varargin)

Expt.Stimvals.et = 'dO';
Expt.Stimvals.e2 = 'ce';
Expt.Stimvals.tf = 0;
Expt.Stimvals.nf = 0;
headonly = 0;
verbose = 0;
j = 1;
while j < nargin
    if strncmpi(varargin{j},'Header',4) %% Read header only
        headonly = 1;
    else
    end
    j = j+1;
end

if verbose tic; end

hname = strrep(name,'.ma','.st');
[ids, vals] = textread(hname,'%s%s');
for j = 1:length(ids)
    if strncmpi(ids(j),'txtfile',5)
        Expt.txtable = strrep(name,'.ma','.txt');
    else
        x = str2num(vals{j});
        if isempty(x) % was a string
            eval(['Expt.Stimvals.' ids{j} ' = vals{j};']);
        else
            eval(['Expt.Stimvals.' ids{j} ' = x;']);
        end
    end
end

if ~isfield(Expt.Stimvals,'fz')
    Expt.Stimvals.fz = 96.04;
end
if ~isfield(Expt.Stimvals,'ce')
    Expt.Stimvals.ce = 1;
end


if headonly
    return;
end

if isfield(Expt.Stimvals,'pt') & Expt.Stimvals.pt > 1
    runseq = 1;
    Expt.isrc = 1;
else
    Expt.isrc = 0;
    runseq = 0;
end

z = textread(name);
if verbose toc; end
ns = 1;
bin = 0;
binw = Expt.Stimvals.bw * 10; %% Timestamp units
Expt.Trials = [];
Expt.Trials(1).Spikes = [];
Expt.Trials(1).Start = 0;
Expt.Trials(1).End = Expt.Stimvals.du;
Expt.Trials(1).Trial = 1;
if Expt.isrc
    Expt.Trials.TrialStart = 0;
end

Expt.Names.dO = 'Orthogonal Disparity';
Expt.Names.ce = 'Correlation';
Expt.Names.or = 'Orientation';
Expt.Names.sf = 'SF';
Expt.Names.tf = 'TF';
Expt.Names.me = 'Eye';
Expt.Names.ip = 'Phase';
Expt.Names.dp = 'dPhase';
Expt.Names.dq = 'dPhase2';
Expt.Names.dx = 'Disparity (H)';
waspreperiod = 0;

Expt.Header.Name = name;
Expt.Header.nsines = 18;
msoff = 0;
bin = 0;
npt = 0;
cols = size(z,2);
for j = 1:size(z,1)
    if runseq & ~isnan(z(j,2)) & ((z(j,4) > 0)  | waspreperiod)
        waspreperiod = 0;
        npt = npt+1;
        Expt.Trials(ns).Start(npt,1) =  Expt.Trials(ns).TrialStart + bin * binw;
        Expt.Trials(ns).id = z(j,4);
        if npt > 1
            Expt.Trials(ns).End(npt,1) =  Expt.Trials(ns).Start(npt,1) + mean(diff(Expt.Trials(ns).Start));
        end
    end
%
% inf inf = new stimulus
    if isinf(z(j,1)) & isinf(z(j,2))
        if bin > 0
            ns = ns+1;
        end
        if runseq
            Expt.Trials(ns).TrialStart = j * 10000;
        end
        if z(j,3) > 0
            Expt.Trials(ns).Start = z(j,3);
        else
            Expt.Trials(ns).Start = j * 10000;
        end
        Expt.Trials(ns).OptionCode = '';
        Expt.Trials(ns).End = Expt.Trials(ns).Start + Expt.Stimvals.du;
        Expt.Trials(ns).Spikes = [];
        Expt.Trials(ns).Trial = ns;
        bin = 0;
        npt = 0;
        Expt.Trials(ns).Nf = Expt.Stimvals.nf;
        if(z(j,4) > 1000)
            Expt.Trials(ns).se = z(j,4);
            Expt.Trials(ns).ls = z(j,4) +2*(Expt.Stimvals.nf-1);
        end
    elseif isnan(z(j,1)) & isnan(z(j,2)) %RespDir
      Expt.Trials(ns).RespDir = z(j,3);
    elseif isinf(z(j,1)) & isnan(z(j,2)) %prestimulus
        for k = 1:z(j,3)
            Expt.Trials(ns).Spikes = [Expt.Trials(ns).Spikes ((bin-1) * binw) -500];
            msoff = msoff+10;
            if msoff > binw
                msoff = 0;
            end
        end
        waspreperiod = 1;
        if z(j,4) 
            Expt.Trials(ns).id = z(j,4);
        end
    elseif isnan(z(j,1)) & isinf(z(j,2)) %postperiod
        test = 1;
    elseif isinf(z(j,1))
        if runseq
% Fill e2 before e1, in case expts are added or e2 is none.
            Expt.Trials(ns).(Expt.Stimvals.e2)(npt,1) = z(j,2);
            Expt.Trials(ns).(Expt.Stimvals.et)(npt,1) = z(j,1);
            Expt.Trials(ns).dO(npt,1) = 0;
            Expt.Trials(ns).ce(npt,1) = 0;
        else
% Inf followed by a number is uncorrelated. Set et and e2 vals, but then
% make sure ce is zero.
            Expt.Trials(ns).(Expt.Stimvals.e2) = z(j,2);
            Expt.Trials(ns).(Expt.Stimvals.et) = 0;
            Expt.Trials(ns).ce = 0;
            Expt.Trials(ns).dO = 0;
            if isfield(Expt.Trials,'me') & isempty(Expt.Trials(ns).me)
                Expt.Trials(ns).me = 0;
            end
        end
        if z(j,2) == -1001
            Expt.Trials(ns).me = -1;
        elseif z(j,2) == -1002
            Expt.Trials(ns).me = 1;
        end
%give spikes timestamps that are within the bin, but semi-randomly
%distributed, by incrementing a 1ms counter each time a spike is
%encountered       
        for k = 1:z(j,3)
            Expt.Trials(ns).Spikes = [Expt.Trials(ns).Spikes msoff+((bin-1) * binw)];
            msoff = msoff+10;
            if msoff > binw
                msoff = 0;
            end
        end
        bin = bin+1;
    else
        if cols > 3 & (z(j,4) | waspreperiod)
            waspreperiod = 0;
            if runseq
                Expt.Trials(ns).ce(npt,1) = Expt.Stimvals.ce;
                Expt.Trials(ns).(Expt.Stimvals.et)(npt,1) = z(j,1);
                Expt.Trials(ns).(Expt.Stimvals.e2)(npt,1) = z(j,2);
                if z(j,1) == -1009 | z(j,2) == -1009
                    Expt.Trials(ns).st(npt,1) = 0;
                else
                    Expt.Trials(ns).st(npt,1) = Expt.Stimvals.st;
                end
                if ~strcmp('me',Expt.Stimvals.et)
                    Expt.Trials(ns).me(npt,1) = 0;
                end
                if z(j,1) == -1005 | z(j,2) == -1005
                    Expt.Trials(ns).ce(npt,1) = 0;
                end
                if z(j,1) == -1002 | z(j,2) == -1002
                    Expt.Trials(ns).me(npt,1) = -1;
                end
                if z(j,1) == -1001 | z(j,2) == -1001
                    Expt.Trials(ns).me(npt,1) = 1;
                end
              
                Expt.Trials(ns).id = z(j,4);
            else
          
                if z(j,1) == -1001
                    Expt.Trials(ns).me = -1;
                elseif z(j,1) == -1002
                    Expt.Trials(ns).me = 1;
                else
                    Expt.Trials(ns).me = Expt.Stimvals.me;
                end
                if z(j,1) == -1005 | z(j,2) == -1005
                    Expt.Trials(ns).st = 0;
                else
                    Expt.Trials(ns).st = Expt.Stimvals.st;
                end
                Expt.Trials(ns).(Expt.Stimvals.e2) = z(j,2);
                if z(j,1) > -999
                    Expt.Trials(ns).(Expt.Stimvals.et) = z(j,1);
                else
                    Expt.Trials(ns).(Expt.Stimvals.et) = 0;
                end

            end
        end
        if (z(j,1) == -1005 | z(j,2) == -1005) & ~runseq
            Expt.Trials(ns).st = 0;
        end
        nspk = z(j,3);
        for k = 1:nspk
            Expt.Trials(ns).Spikes = [Expt.Trials(ns).Spikes msoff + ((bin-1) * binw)]; 
            msoff = msoff+10;
            if msoff > binw
                msoff = 0;
            end
        end
        bin = bin+1;
   end
    
end

for j = 1:length(Expt.Trials)
    if(isfield(Expt.Trials,'st') & isempty(Expt.Trials(j).st))
        Expt.Trials(j).st = Expt.Stimvals.st;
    end
    if(isfield(Expt.Trials,'ce') & isempty(Expt.Trials(j).ce))
        Expt.Trials(j).ce = Expt.Stimvals.ce;
    end
    Expt.Trials(j).Count = sum(Expt.Trials(j).Spikes > 0);
   Expt.Trials(j).Spikes = Expt.Trials(j).Spikes';  %% Make a colimn, like offline
end
if verbose toc; end