function TestNev2Mat(datafolder, savefolder)

foldernames = ls(datafolder);
foldercounts = length(foldernames);

for i=3:foldercounts-1
    clear NEV;
    foldername = foldernames(i,1:end-1);
    fullpath = [datafolder '\' foldername '\' foldername 'a.nev'];
    NEV = openNEV('read', 'report', fullpath);
    
    nfullpath = [savefolder '\' foldername '.mat'];

    save(nfullpath, 'NEV');
end