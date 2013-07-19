%% Convert a NEV file into a MAT file for fast data read and smaller file
%  size
%
%  Use NEV2MAT

function NEV2MAT

datafolder = uigetdir('D:\', 'Select the folder containing Cerebus folders...');
savefolder = uigetdir('D:\', 'Select the folder where you would like to save MAT files...');
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