function [probit, blocklist]= readpsychsum(file, varargin)

j = 1;
k = 1;
doall = 0;
skipblock = 0;
lastblock = 10000000;
collapsedata = 0;
blockstouse = [];
useblocks = [];
nc = 0;
showtimes = 1;
NOTSET = -100000;

varargon = {};
conditions = {};
while(j < nargin)
    if(strncmpi(varargin{j},'all',3))
        doall = 1;
    elseif(strncmpi(varargin{j},'collapse',4))
      collapsedata = 1;
    elseif(strncmpi(varargin{j},'skip',4))
	j = j+1;
        skipblock = varargin{j};
    elseif(strncmpi(varargin{j},'blocks',4))
        j = j+1;
   %blockstouse is a vector of 0s or 1s describing which blocks to use.
        blockstouse = varargin{j};
    elseif(strncmpi(varargin{j},'useblocks',4))
        j = j+1;
        useblocks = varargin{j};
    elseif(strncmpi(varargin{j},'sort',4))
	j = j+1;
    nc = nc+1;
    conditions{nc}  = varargin{j};
    clabel{nc} = strrep(varargin{j},'=','');
    nvals(nc) = 0;
    elseif(strncmpi(varargin{j},'last',4))
	j = j+1;
        lastblock = varargin{j};
    else
        varargon{k} = varargin{j};
        k = k+1;
    end
    j = j+1;
end

expno = 1;
probit = [];
ndat = 1;
gotdat = 0;
vals = [];
allvals = [];
valset = [];
fid = fopen(file,'r');
tline = fgetl(fid);
nblocks = 0;
nline = 1;
if showtimes
tic;
end
while(ischar(tline))
  [x , n]  = sscanf(tline,'%f %d %d');
  if(n > 2)
    if(x(2) > 0)
      probit(ndat).x = x(1);
      probit(ndat).n = x(2);
      probit(ndat).p = x(3)/x(2);
      probit(ndat).resp = x(3);
      probit(ndat).expno = expno;
      probit(ndat).name = name;
      probit(ndat).block = nblocks;
      probit(ndat).i = ndat;
      for j = 1:nc
          if j > length(valset) | ~valset(j)
              vals(j) = NOTSET;
          end
          eval(['probit(ndat).' clabel{j} '= vals(j);']);
      end
      probit(ndat).vals = vals;
      blocks(ndat) = nblocks;
      ndat = ndat+1;
      gotdat = gotdat +1;
    end
  end
  if(strncmpi(tline,'Experiment',10))
    expno = sscanf(tline,'Experiment %d');
    j = strfind(tline,':');
    if ~isempty(j)
      name = tline(j(1)+1:end);
    end
    gotdat = 0;
      nblocks = nblocks +1;
      blocks(ndat) = nblocks;
      blockn(ndat) = ndat;
      blockstart(nblocks) = ndat;
      blocknames{ndat} = name;
  elseif(strncmpi(tline,'Run ended',8))
      blocktimes{blockstart(nblocks)} = tline(10:end);
      if length(vals) > 0
          allvals(:,nblocks) = vals;
      end
      valset = [];
  elseif(strncmpi(tline,'Condition',8) | strncmpi(tline,'Step',4) | strncmpi(tline,'Stimulus',6))
    if(gotdat & ~collapsedata)
      expno = expno+100;
    end
    j = strfind(tline,':');
    if ~isempty(j)
      if(gotdat)
	name = tline(j(1)+1:end);
      else
	name = [name tline(j(1)+1:end)];
      end
    end
    for j = 1:nc
        id = strfind(tline,[conditions{j} '=']);
        if ~isempty(id)
            vals(j) = sscanf(tline(id+length(conditions{j})+1:end),'%f');
            valset(j) = 1; 
            nvals(j) = nvals(j)+1;
%            allvals(j,nvals(j)) = vals(j);
        end
    end
    
  end
  tline = fgetl(fid);
  nline = nline+1;
end
fclose(fid);
if showtimes
toc;
end

id = length(blockn)+1:length(probit);
for j = id
    blockn(j) = 0;
    blocknames{j} = 0;
end
id = find(blockn > 0);
if id(end) > length(probit)
  id = id(1:end-1);
end

blocks = blocks(1:length(probit));
blocklist.blocks = blocks(id);
blocklist.blockids = id;
blocklist.names = {blocknames{id}};
blocklist.n = [probit(id).n];
blocklist.times = {blocktimes{id}};
blocklist.alldata = probit;



idall = id(1:end);
idall(end+1) = length(probit);

for j = 1:(length(id));
    prange = idall(j):(idall(j+1)-1);
  blocklist.n(j) = sum([probit(prange).n]);
  for k = 1:length(conditions)
      eval(['blocklist.' conditions{k} '{j} = unique([probit(prange).' conditions{k} ']);']);
  end
end



if showtimes
toc;
end

ndat = 1;
allprobit = [];

if isempty(useblocks)
if isempty(blockstouse)
  if lastblock < 1
    useblocks = skipblock:lastblock;
  else
    useblocks = skipblock:max(blocks);
  end
else
  useblocks = blocklist.blocks(find(blockstouse > 0));
end
end

for j = unique([probit.expno])
  idx = find([probit.expno] == j);
  for val = unique([probit(idx).x])
    id = find([probit.x] == val & [probit.expno] == j & ismember(blocks,useblocks));
    if ~isempty(id)
      allprobit(ndat).x = val;
      allprobit(ndat).expno = j;
      allprobit(ndat).n = sum([probit(id).n]);
      allprobit(ndat).resp = sum([probit(id).resp]);
      allprobit(ndat).p = allprobit(ndat).resp/allprobit(ndat).n;
      allprobit(ndat).name = probit(id(end)).name;
      ndat = ndat+1;
    end
  end
end


if nc > 0
for j =1:length(probit)
    probitvals(:,j) = probit(j).vals';
end
end
ndat = 1;
if size(allvals,1) == 1
    j = 1;
    for sval = unique(allvals(j,:))
        idx = find([probit.vals] == sval);
        for val = unique([probit(idx).x])
            id = find([probit.x] == val & [probit.vals] == sval & ismember(blocks,useblocks));
            if ~isempty(id)
                sortprobit(ndat).x = val;
                sortprobit(ndat).expno = j;
                sortprobit(ndat).n = sum([probit(id).n]);
                sortprobit(ndat).resp = sum([probit(id).resp]);
                sortprobit(ndat).p = sortprobit(ndat).resp/sortprobit(ndat).n;
                sortprobit(ndat).name = probit(id(end)).name;
                sortprobit(ndat).vals = sval;
                eval(['sortprobit(ndat).' clabel{1} '= sval;']);
                ndat = ndat+1;
            end
        end
        j = j+1;
    end
end


if size(allvals,1) > 1

    j = 1;
    ncon = 1;
    for dim = 1:size(allvals,1)
        chkvals{dim} = unique(allvals(dim,:));
        ncon = ncon * length(chkvals{dim});
    end
    for j = 1:ncon
        dm = 1;
        for dim = 1:size(allvals,1)
            k = mod(floor((j-1)/dm),length(chkvals{dim}))+1;
            idim{dim} = find(probitvals(dim,:) == chkvals{dim}(k));
            dm = dm * length(chkvals{dim});
            if(dim == 1)
                idu = idim{1};
            else
                idu = intersect(idu,idim{dim})
            end
            chklist(j,dim) = k;
        end
        idlist{j} = idu;
        idx = idu;
        tmp = probit(idx);
        for val = unique([probit(idx).x])
            id = find([probit.x] == val & ismember([probit.i],idx) & ismember(blocks,useblocks));
            if ~isempty(id)
                sortprobit(ndat).x = val;
                sortprobit(ndat).expno = j;
                sortprobit(ndat).n = sum([probit(id).n]);
                sortprobit(ndat).resp = sum([probit(id).resp]);
                sortprobit(ndat).p = sortprobit(ndat).resp/sortprobit(ndat).n;
                sortprobit(ndat).name = probit(id(end)).name;
%%                sortprobit(ndat).vals = [aval bval];
                sortprobit(ndat).id = id;
%%                eval(['sortprobit(ndat).' clabel{1} '= aval;']);
  %%              eval(['sortprobit(ndat).' clabel{2} '= bval;']);
                sortprobit(ndat).sname = sprintf('%d:',j);

                for k = 1:length(clabel);
                    eval(['sortprobit(ndat).' clabel{k} '= mean([probit(id).' clabel{k} ']);']);
                    a = eval(['sortprobit(ndat).' clabel{k}]);
                    sortprobit(ndat).sname = [sortprobit(ndat).sname sprintf(' %s=%f',clabel{k},a)];
                end
                
                ndat = ndat+1;
            end
        end
    end
    if 0
    for aval = unique(allvals(1,:))
        for bval = unique(allvals(2,:))
        idx = find(probitvals(1,:) == aval & probitvals(2,:) == bval);
        tmp = probit(idx);
        for val = unique([probit(idx).x])
            id = find([probit.x] == val & probitvals(1,:) == aval & probitvals(2,:) == bval & ismember(blocks,useblocks));
            if ~isempty(id)
                sortprobit(ndat).x = val;
                sortprobit(ndat).expno = j;
                sortprobit(ndat).n = sum([probit(id).n]);
                sortprobit(ndat).resp = sum([probit(id).resp]);
                sortprobit(ndat).p = sortprobit(ndat).resp/sortprobit(ndat).n;
                sortprobit(ndat).name = probit(id(end)).name;
                sortprobit(ndat).vals = [aval bval];
                sortprobit(ndat).id = id;
                eval(['sortprobit(ndat).' clabel{1} '= aval;']);
                eval(['sortprobit(ndat).' clabel{2} '= bval;']);
                sortprobit(ndat).sname = sprintf('%s=%f %s=%f',clabel{1},aval,clabel{2},bval);
                for k = 3:length(clabel);
                    eval(['sortprobit(ndat).' clabel{k} '= mean([probit(id).' clabel{k} ']);']);
                    a = eval(['sortprobit(ndat).' clabel{k}]);
                    sortprobit(ndat).sname = [sortprobit(ndat).sname sprintf(' %s=%f',clabel{k},a)];
                end
                
                ndat = ndat+1;
            end
        end
        j = j+1;
    end
    end
end

end

if(nc)
    probit = sortprobit;
else
    probit = allprobit;
end

toc;