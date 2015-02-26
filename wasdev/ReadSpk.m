function SPK = ReadSpk(spkfile);


mafile = strrep(spkfile,'.spk','.mat');


cluster = (max(find(spkfile=='_')));
if cluster
  cluster = str2num(spkfile(cluster+1));
  spkfile = strrep(spkfile,sprintf('_%d',cluster),'');
else
  cluster = 1;
end

if (~exist(spkfile))
  fprintf('No File%s',spkfile);
  return;
end

    fprintf('Reading %s\n',spkfile);
  infid = fopen(spkfile,'r');
  spiketimes = [];
  SPK.spikes = {};
  SPK.strings = {};
  stimcount = 0; % gets incremented first
  while ~feof(infid)
    scrap = fgets(infid);
    if(scrap ~= -1)
      clusterno = sscanf(scrap,'%f');
      if ~isempty(clusterno) %then the line begins with a digit
	ldata = sscanf(scrap,'%f');
	if(ldata(1) == cluster)
	  spiketimes = [spiketimes ldata(2)];
	  SPK.Trials(stimcount).Spikes = spiketimes;
	end
      elseif(strncmp(scrap,'End',3))
	stimcount = stimcount +1;
      elseif(strncmp(scrap,'Stim',4))
	spiketimes = [];
	if(stimcount == 0 | stimcount <= length(SPK.Trials))
	  stimcount = stimcount +1;
	end
	SPK.Trials(stimcount).string = scrap;
	[SPK.Trials(stimcount).Start] = sscanf(scrap,'Stim %d');
	dur = sscanf(scrap,'Stim %*d, %f');
	SPK.Trials(stimcount).End = SPK.Trials(stimcount).Start + ...
	    10000 * dur;
	idx = findstr(scrap,'sf=');
	if(length(idx))
	  SPK.Trials(stimcount).sf  = sscanf(scrap(idx(1):end),'sf=%f');
	end
      end
    end
  end
  SPK.ntrials = length(SPK.strings);
  SPK.name = mafile;
  save(mafile,'SPK');

