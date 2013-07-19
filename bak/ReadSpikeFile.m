function Spikes = ReadSpikeFile(spkfile, varargin)

loadallspikes = 0;
j = 1;
while j  <= length(varargin)
    if strncmpi(varargin{j},'allprobes',6)
        loadallspikes = 1;
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
load(spkfile);
if size(Spikes.values,2) > 100
    Spikes.values = Spikes.values';
end
    Spikes.times = Spikes.times .* 10000;
    Spikes.times = reshape(Spikes.times,length(Spikes.times),1);
    if size(Spikes.codes,2) == 1
        Spikes.codes(:,2) = Spikes.codes(:,1);
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
if loadallspikes
    xname = regexprep(spkfile,'.p([0-9])*t','.p$1xt');
    if exist(xname)
    X = load(xname);
    Spikes.xmaxv = X.Spikes.maxv;
    Spikes.xvalues = X.Spikes.values;
    Spikes.xchans = X.Spikes.chspk;
    end
end