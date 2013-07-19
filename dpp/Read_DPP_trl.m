fid = fopen('/bgc/bgc/data/dufus/299/duf299.0.2grating.DPP.trl','r');
if fid==-1
  return
end

stimcount = 1;
while ~feof(fid)
    scrap = fgets(fid);
    trialno = sscanf(scrap,'%f');
    if ~isempty(trialno) %then the line begins with a digit
	trialvals = sscanf(scrap, '%d(%d,%f): %*s %*s %d %f');
	stimdescript = sscanf(scrap,'%*s %*s %s',1);
            while length(stimdescript)>0
                [param,rem]=strtok(stimdescript,'=');
                val = sscanf(rem(2:end),'%f');
                % Work out what param this value relates to and record it:
                switch strrep(param,',','') %removeunwanted commas
		 case 'lxo'
		  TC.grating.eye = 0;
		 case 'rxo'
		  TC.grating.eye = 1;
		  eye = R;
		 case 'bxo'
		  TC.grating.eye = 2;
		 case 'sf'
		  TC.grating.orientation(stimcount) = val;
		 case 'dp'
		  TC.grating.phasedisp(stimcount) = val;
		 case 'dq'
		  TC.grating.phasebdisp(stimcount) = val;
		end
		                [rem,stimdescript] = strtok(rem,',');
            end % finish going through stim description
	    stimcount = stimcount+1;
    end
end


