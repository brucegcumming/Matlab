function FullV = LoadFullV(name, varargin)
%FullV = LoadFullV(name)
%loads a FullV file from disk, and performs any necessary scaling
%

if iscell(name)
    for j = 1:length(name)
        res(j) = LoadFullV(name{j}, varargin{:});
    end
    return;
end

highpass = 1;
smoothw = 100;
keepsmooth = 0;
convert = 1;
toint = 0;

FullV = [];
if ~exist(name,'file')
    fprintf('%s Does not exist\n',name);
    return;
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'highpass',5)
        j = j+1;
        smoothw = varargin{j};
        keepsmooth = 1;
    elseif strncmpi(varargin{j},'checkint',6)
        toint = 3;
    elseif strncmpi(varargin{j},'converttoint',5)
        toint = 1;
    elseif strncmpi(varargin{j},'saveint',5)
        toint = 2;
    elseif strncmpi(varargin{j},'noconvert',5)
        convert = 0;
    end
    j = j+1;
end


isint = 0;
ts = now;
load(name);
FullV.loadname = name;
if toint && isfloat(FullV.V)
    if isinteger(FullV.V) %alreay i sint
        isint = 1;
        fprintf('%s isalready ints\n',name);
        return;
    elseif toint == 3 %just checking
        fprintf('%s is double\n',name);
        return;
    end
    intscale(2) = 32000;
    intscale(1) = max(abs(FullV.V(:)));
     FullV.V = int16(round(double(FullV.V .* intscale(2)./intscale(1))));    
     FullV.intscale = intscale;
     if toint == 2
         fprintf('Saving %s as int\n',name);
         save(name,'FullV');
     end
    return;
end

method = 3;  %change manually for testing
FullV.initialloadtime = mytoc(ts);
d = whos('FullV');
FullV.readrate = d.bytes./FullV.initialloadtime;
FullV.size = d.bytes./(1024.*1024);
if convert && isinteger(FullV.V)  && isfield(FullV,'intscale')
    if method == 1  %try to save memory
        np = size(FullV.V,1);
        for j = 1:np
            NewV(np+1-j,:) = double(FullV.V(np+1-j,:)) .* FullV.intscale(1)/FullV.intscale(2);
            FullV.V = FullV.V(1:np-j,:);
        end
        FullV.V = NewV;
    elseif method == 2
        FullV.V = double(FullV.V) .*FullV.intscale(1)/FullV.intscale(2);
    elseif method == 3
        FullV.V = double(FullV.V);
        np = size(FullV.V,1);
        for j = 1:np
            FullV.V(j,:) = FullV.V(j,:) .* FullV.intscale(1)/FullV.intscale(2);
        end
    else
        NewV = double(FullV.V) .*FullV.intscale(1)/FullV.intscale(2);
        FullV.V = NewV;
        clear NewV;
    end
end
FullV.loadtime = mytoc(ts);
if keepsmooth
    FullV.highpass = smoothw;
end
if isfield(FullV,'chspk') &&  ~isempty(strfind(name,'p1FullV')) && FullV.chspk == 96
    FullV.chspk = 1;
end
if isfield(FullV,'highpass') && ~isnan(FullV.highpass)
    ts = now;
    if FullV.highpass > 0
        smoothw = FullV.highpass;
    end
     sm = smooth(double(FullV.V),smoothw);
     if keepsmooth
         FullV.Vsmooth = sm;
     end
     FullV.V = int16(round(double(FullV.V) - sm));
     FullV.filtertime = mytoc(ts);
end
