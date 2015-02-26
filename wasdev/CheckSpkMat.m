function [result, details] = CheckSpkMat(varargin)

% CheckSpkMat(varargin)  checks contents of matlab variables from Spike3
% matlab files.
% CheckSpkMat(name,.....) loads file named first
% CheckSpkMat(Chx, Chy ....) Checks the matlab variables that are passed.
checkmode = 1;
result = 1;
details = [];
filename = [];
j = 1;
while j <= length(varargin)
    v = varargin{j};
    if isstruct(v)
        if isfield(v,'title')
            if strfind(v.title,'Spike')
                ChSpk = v;
            elseif strfind(v.title,'StimOn')
                ChStim = v;
            elseif strfind(v.title,'Serial')
                ChSer = v;
            elseif strfind(v.title,'Mains')
                ChMains = v;
            elseif strfind(v.title,'VTR')
                ChVTR = v;
            end
        end
    elseif ischar(varargin{j}) & exist(varargin{j},'file')
        filename = varargin{j};
        A = load(varargin{j});
        f = fields(A);
        for k = 1:length(f)
            if isfield(A.(f{k}),'title')
                if strfind(A.(f{k}).title,'Mains')
                    ChMains = A.(f{k});
                end
            end
        end
    elseif strncmpi(varargin{j},'mains',5)
        checkmode = 2;
    elseif strncmpi(varargin{j},'vtr',3)
        checkmode = 3;
    end
    j = j+1;
end
 
if checkmode == 1
    bstimes = ChStim.times(find(ChStim.level == 1));
    estimes = ChStim.times(find(ChStim.level == 0));
    sid = strmatch('bss',ChSer.text);
    stimes= ChSer.times(sid);
    eid = strmatch('ess',ChSer.text);
    etimes= ChSer.times(eid);
    for j = 2:length(etimes)-1;
        id = find(estimes < etimes(j)+0.015);
        delays(j,3) = etimes(j)-estimes(id(end));
        delays(j,4) = etimes(j)-estimes(id(end)+1);

        id = find(stimes < etimes(j));
        k = id(end);
        id = find(bstimes < stimes(k)+0.015);
        delays(j,1) = stimes(k)-bstimes(id(end));
        
        delays(j,2) = stimes(k)-bstimes(id(end)+1);
        starts(j) = stimes(k);
    end
    plot(delays(:,1),delays(:,2),'o')
    details.delays = delays;   
    details.starts = starts;
elseif checkmode == 2 & exist('ChMains')
    [a, b] = min(diff(ChMains.times));
    if a < 0.01
        fprintf('%s has short mains interval (%.6f) at %.2f\n',filename,a,ChMains.times(b));
        result = 0;
        id = find(diff(ChMains.times) < 0.01);
        details.times = ChMains.times(id);
    end
    details.minmains = a;
elseif checkmode == 3 & exist('ChVTR')
    [a, b] = min(diff(ChVTR.times));
    if a < 0.005
        fprintf('%s has short VTR interval (%.6f) at %.2f\n',filename,a,ChVTR.times(b));
        result = 0;
        id = find(diff(ChVTR.times) < 0.01);
        details.times = ChVTR.times(id);
    end
    details.minvtr = a;
else
    result = 1;
end

