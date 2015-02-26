function simplesave(somefilename, laps)

deleteprev = 0; %without this causes trouble wiht Maveriks and Samba
fprintf('Testing %s\n',somefilename);
backupname = strrep(somefilename,'.mat','bak.mat');
for i = [1:laps]
    kashk.r = randi(10, [10000000,1]); 
%either method of deleting the file first seems to prevent the problems
    if exist(somefilename) && deleteprev == 1
        delete(somefilename);
    elseif exist(somefilename) && deleteprev == 2
        movefile(somefilename,backupname);
    end
    fprintf('Saving %d: %s at %s\n',i,somefilename,datestr(now));
%-v7.3  produces error at write time.
%-v6  gets rid of error
    save(somefilename,  'kashk');
    fprintf('Loading %s at %s\n',somefilename,datestr(now));
    load(somefilename); 
end

