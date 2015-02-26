function CheckLists(DATA)

pairs = {'ORBW' '/bgc/bgc/anal/orbw/lemORBW.lst' ; ...
'OTRC' '/bgc/bgc/anal/orbw/otrc.lst'};

[a,dp] = splitpath(DATA.datafilename);
d = dir(dp);
nf = 0;
for j = 1:length(d)
if strfind(d(j).name,'.mat') 
for k = 1:size(pairs,1)
if strfind(d(j).name,pairs{k,1})
if isempty(cmb.CheckListForName(pairs{k,2},d(j).name))
nf = nf+1;
dat.files{nf} = [dp '/' d(j).name];
dat.lists{nf} = pairs{k,2};
else
fprintf('%s Already in %s\n',d(j).name,pairs{k,2});
end
dat.nf = nf;
end
end
end
end
if nf
SPACE = 10;
cw=8;
ch=10;
figure('Position',[10 10 cw*60, (ch+SPACE) * (nf+2)]);
fp = get(gcf,'Position');
bp = [10 10 fp(4) 20];
for j = 1:nf
bp = [10 bp(2)+bp(4)+SPACE fp(3) ch * 2];
uicontrol(gcf,'Style', 'checkbox',...
'String', [dat.files{j} dat.lists{j}], 'Tag', 'SuffList', 'Position', bp,'UserData',j);
end
bp = [10 10 cw * 5 20];
uicontrol(gcf,'Style', 'pushbutton',...
'String','go', 'callback',@cmb.cplists);
set(gcf,'UserData',dat);

end


