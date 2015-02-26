function ok = CheckSpikeFile(C, fullvname)
%Check spk files are up to date for a cluster
    ok = 1;
    if ~isfield(C,'spkfile')
        fprintf('No spike file in Cluster %d\n',C.probe(1));
        return;
    end
    [monkey, a,b]  = GetMonkeyName(C.spkfile);
    id = findstr(fullvname,monkey);
    bid = findstr(C.spkfile,monkey);
    spkfile = [fullvname(1:id(1)) C.spkfile(1+bid(1):end)];
    spkfile = strrep(spkfile,['/Spikes/' monkey],'/Spikes/');

    ddir = fileparts(fullvname);
    [a,fname] = fileparts(C.spkfile);
    %if spkfile is already new style in cluster, newspkfile is nonsense
    spkfile = [ddir '/Spikes/' fname '.mat'];
    newspkfile = strrep(spkfile,'/Spikes/',['/Spikes/' monkey]);
    if ~exist(spkfile) && ~exist(newspkfile)
        fprintf('Missing Spk file %s and %s\n',spkfile,newspkfile);
        ok = 0;
    else
        if exist(newspkfile)
            spkfile = newspkfile;
        end
        d = dir(spkfile);
        fprintf('%s dir finds %d matches\n',spkfile,length(d));
        savetime = C.savetime(1);
        if d.datenum < savetime-0.01
            ok = 0;
            fprintf('Spk file %s is older (%s) than Cluster (%s)\n',spkfile,d.date,datestr(C.savetime(1)));
        elseif d.datenum > savetime +1
%            fprintf('Spk file %s is much newer (%s) than Cluster (%s)\n',spkfile,d.date,datestr(C.savetime(1)));
            load(spkfile);
            if size(Spikes.values,1) ~= C.nspks
                ok = 0;
                fprintf('End event list mismatch %d vs %d\n',size(Spikes.values,1),C.nspks);
            end
        end
    end