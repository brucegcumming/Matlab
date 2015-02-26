function mupdate(name)
%update network/local copy of file from local/network copy

src = which(name);
if isempty(regexp(src,'^[A-F]:')) %network file first
    tgt = regexprep(src,'^[Y-Z]','C');
else
    tgt = regexprep(src,'^[A-F]','Z');
end

srcd = dir(src);
tgtd = dir(tgt);
if srcd.datenum < tgtd.datenum
    a = src;
    src = tgt;
    tgt = a;
    a = srcd;
    srcd = tgtd;
    tgtd = a;
    s = sprintf('Copy %s (newer) to %s (Path)\n%s %s\n%s %s\n',src(1:2),tgt(1:2),src,srcd.date,tgt,tgtd.date);
elseif srcd.datenum == tgtd.datenum
    s = sprintf('Copy %s to %s (same age)\n%s %s\n%s %s\n',src(1:2),tgt(1:2),src,srcd.date,tgt,tgtd.date);
else
    s = sprintf('Copy %s (newer, in path) to %s\n%s %s\n%s %s\n',src(1:2),tgt(1:2),src,srcd.date,tgt,tgtd.date);
end
ok = questdlg(s,'Source File Copy','Copy','Cancel', 'Copy');
if strcmp(ok,'Copy');
    fprintf('copying %s to %s\n',src,tgt);
    copyfile(src,tgt);
end
    