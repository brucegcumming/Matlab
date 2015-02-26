function [DATA, Stimvals, needv, retval] = CheckCombine(DATA, interactive)

if nargin < 2
interactive = 1;
end
retval = 1;
err = '';
id = get(DATA.elst,'value');
for j = 1:length(id)
c.stims(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'st');
c.bws(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'ob');
c.mes(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'me');
c.tfs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'tf');
c.sfs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'sf');
c.ors(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'or');
c.sls(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'sl');
c.ces(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'ce');
c.xos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'xo');
c.yos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'yo');
c.fxs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'fx');
c.fys(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'fy');
c.jxs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'jx');
c.cos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'co');
c.ods(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'od');
c.Bcs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'Bc');
c.bos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'bo');
c.backxos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'backxo');
c.backyos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'backyo');
ets{j} = DATA.Expts{DATA.expid(id(j))}.Stimvals.et;
e2s{j} = DATA.Expts{DATA.expid(id(j))}.Stimvals.e2;
Bs = GetEval(DATA.Expts{DATA.expid(id(j))},'Bs');
if ischar(Bs)
a = strmatch(Bs,DATA.stimnames);
if length(a)
c.Bss(j) = DATA.stimnamecodes(a);
else
c.Bss(j) = NaN;
end
else
c.Bss(j) = Bs;
end
Stimvals = DATA.Expts{DATA.expid(id(j))}.Stimvals;
end
%check if ori of back stim matters
bstim = sum(c.Bss ==2 & c.sls > 0) + sum(ismember(c.Bss ,[3 4 11]));



% make c.(xx) a list of unigue values of xx
% if there are > 1 of these, it needs to be a field in Trials
f = fields(c);
needv = {};
n = 1;
for j = 1:length(f)
vals = unique(c.(f{j}));
vals = vals(find(~isnan(vals)));
c.(f{j}) = vals;
if length(vals) > 1
needv{n} = f{j}(1:end-1); %remove the 's'
n = n+1;
end
end

Stimvals.BlockedStim = 0;
for j = 1:length(e2s)
if isempty(e2s{j}) 
e2s{j} = 'e0';
end
if isempty(ets{j}) 
ets{j} = 'e0';
end
end
blank = strmatch('none',c.stims);
nstims = length(c.stims) - length(blank);
blank = strmatch('e0',unique(e2s));
e2lst = unique(e2s);
e1lst = unique(ets);
ne2 = length(e2lst) - length(blank);
bws = unique(c.bws);
bws = bws(find(~isnan(bws)));
xos = unique(c.xos);
yos = unique(c.yos);
ods = unique(c.ods);
%If there are blocks with different values for some paramerter, combine them
%and be sufe to fill trials, then set to expt 2. 
%Careful - this forces combine to generate a new name. 
if length(bws) > 1 & ne2 == 0 & Stimvals.st== 21 %bw only means anything for stim image
err = [err sprintf('%d Bandwidths in %s',length(unique(bws)),DATA.outname) sprintf('%.2f ',unique(bws))];
for j = id(1:end)
DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},'ob');
Stimvals.e2 = 'ob';
Stimvals.BlockedStim = 1;
end
elseif ne2 == 1 & length(blank) > 0 &  strmatch('me',e2lst) %Expt with me varying + one without
for j = id(1:end)
DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},'me');
Stimvals.e2 = 'me';
Stimvals.BlockedStim = 1;
end        
end
if Stimvals.st == 2 && sum(c.sls > 0) == 0 %don't check or for RDS with sl=0
rdssl =1;
else
rdssl = 0;
end

%don't check SF for RDS/RLS/Cylinder
if length(c.sfs) >1 & isempty(strmatch('sf',{e1lst{:} e2lst{:}}))  & ~ismember(Stimvals.st,[2 15 11])
err = [err sprintf('%d SFS in %s',length(c.sfs),DATA.outname) sprintf('%.2f ',c.sfs)];
end
if length(c.ors) > 1 & isempty(strmatch('or',{e1lst{:} e2lst{:}})) & rdssl == 0
err = [err sprintf('%d Oris in %s',length(c.ors),DATA.outname) sprintf('%.2f ',c.ors)];
end
if length(c.ods) > 1 & isempty(strmatch('od',{e1lst{:} e2lst{:}})) & rdssl == 0
err = [err sprintf('%d Ori diffs in %s',length(c.ors),DATA.outname) sprintf('%.2f ',c.ods)];
end
if length(c.sls) > 1 & isempty(strmatch('sl',{e1lst{:} e2lst{:}}))
err = [err sprintf('%d Seed Loops in %s',length(c.sls),DATA.outname) sprintf('%.2f ',c.sls)];
end
if length(c.jxs) > 1 & isempty(strmatch('jx',{e1lst{:} e2lst{:}}))
err = [err sprintf('%d jx vals in %s',length(c.jxs),DATA.outname) sprintf('%.3f ',c.jxs)];
end
if range(xos) > 0.5 & isempty(strmatch('Op',{e1lst{:} e2lst{:}})) & isempty(strmatch('Pp',{e1lst{:} e2lst{:}}))
err = [err sprintf('%d X pos (%.1f deg) in %s',length(c.xos),range(xos),DATA.outname) sprintf('%.2f ',c.xos)];
end

if bstim > 1 && length(c.bos) > 1 && isempty(strmatch('b',{e1lst{:} e2lst{:}}))
err = [err sprintf('%d BackOr in %s',length(c.bos),DATA.outname) sprintf('%.3f ',c.bos)];
end
f = {'cos' 'backxos' 'backyos' 'Bcs'};
ex = {'co' 'backxo' 'backyo' 'Bc'};
names = {'Contrast' 'Back Xpos' 'Back Ypos' 'Back Contrast'};
for j = 1:length(f)
if length(c.(f{j})) > 1 & isempty(strmatch(ex{j},{e1lst{:} e2lst{:}}))
err = [err sprintf('%d %s in %s',length(c.(f{j})),names{j},DATA.outname) sprintf('%.3f ',c.(f{j}))];
end
end
if length(err) > 1
if interactive 
a = questdlg(err,'Combine these Files','Cancel','OK','OK');
if strmatch(a,'Cancel')
retval = -1;
end
else
DATA = AddError(DATA,err);
end
end


