function Spikes = ReadSpikeFile(spkfile, varargin)
%Spikes = ReadSpikeFile(spkfile, varargin)
%    ....,'allprobes')  loads in other probes from xspk file
%    ....,'double')  converts int to double 

loadallspikes = 0;
converttodouble = 0;
j = 1;
while j  <= length(varargin)
    if strncmpi(varargin{j},'allprobes',6)
        loadallspikes = 1;
    elseif strncmpi(varargin{j},'double',6)
        converttodouble = 1;
    end
    j = j+1;
end
spkpath = ['Spikes/' GetMonkeyName(spkfile)];



if ~exist(spkfile)
    spkfile = strrep(spkfile,['Spikes/' GetMonkeyName(spkfile)],'Spikes/');
    if ~exist(spkfile)
        fprintf('Cannot read %s\n',spkfile);
        Spikes = [];
        return;
    end
end
ts = now;
a = load(spkfile);
if isfield(a,'Spikes')
    Spikes = a.Spikes;
else
    f = fields(a);
    if length(f) ==1
        Spikes = a.(f{1});
    end    
end
if size(Spikes.values,2) > 100
    Spikes.values = Spikes.values';
end
    Spikes.times = Spikes.times .* 10000;
    Spikes.times = reshape(Spikes.times,length(Spikes.times),1);
    if size(Spikes.codes,2) == 1
        Spikes.codes(:,2) = Spikes.codes(:,1);
    end
    Spikes.Header.bysuffix = 0;
    if ~isfield(Spikes.Header,'matfile')
        fprintf('SpkFile %s missing Header.matfile\n',spkfile);
        Spikes.Header.matfile = '';
    end
    if regexp(Spikes.Header.matfile,'\.[0-9]*\.mat')
        Spikes.Header.bysuffix = 1;
    end
        
    if ~isfield(Spikes,'maxint')
        Spikes.maxint = 32000;
        if strmatch(class(Spikes.values),'int16')
            if ~isfield(Spikes,'maxv')
                fprintf('Integer Spikes without maxv!!\n');
                Spikes.maxv = double(max(abs(Spikes.values(:)))./5);
            end
            Spikes.VRange = double([min(Spikes.values(:))  max(Spikes.values(:))]).*Spikes.maxv./Spikes.maxint;
        else
            Spikes.maxint = 1;
            Spikes.maxv = max(abs(Spikes.values(:)));
            Spikes.VRange = [min(Spikes.values(:))  max(Spikes.values(:))];
        end
    else
        Spikes.VRange = double([min(Spikes.values(:))  max(Spikes.values(:))]).*Spikes.maxv./Spikes.maxint;
    end
    Spikes.Vscale = Spikes.maxv./Spikes.maxint;
    if isfield(Spikes,'Vrange') && isempty(Spikes.Vrange)
        
    end
if converttodouble && isinteger(Spikes.values)
    Spikes.values = double(Spikes.values) .* Spikes.Vscale;
end
if loadallspikes
    xname = regexprep(spkfile,'.p([0-9])*t','.p$1xt');
    if exist(xname)
        X = load(xname);
        Spikes.xmaxv = X.Spikes.maxv;
        if ndims(X.Spikes.values) == 2
            Spikes.xvalues(1,:,:) = X.Spikes.values;
        else
            Spikes.xvalues = X.Spikes.values;
        end
        if ~isfield(X.Spikes,'xVrange')
            Spikes.xVrange = double([min(X.Spikes.values(:)) max(X.Spikes.values(:))]) .* X.Spikes.maxv./X.Spikes.maxint;
        else
            Spikes.xVrange = X.Spikes.xVrange;
        end
        Spikes.xVscale = X.Spikes.maxv./X.Spikes.maxint;
        Spikes.xchans = X.Spikes.chspk;
        if isfield(X.Spikes,'TriggerV')
            Spikes.TriggerV = X.Spikes.TriggerV;
        end
        if converttodouble && isinteger(Spikes.xvalues)
            Spikes.xvalues = double(Spikes.xvalues) .* Spikes.xVscale;
        end

    end
end
Spikes.Header.loaddur = mytoc(ts);
Spikes.Header.loadtime = ts;
