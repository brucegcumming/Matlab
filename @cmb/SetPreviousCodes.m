function [DATA, done] = SetPreviousCodes(DATA, varargin)
%if online and have stored codes for this probe/expt, put them in spks
done = 0;
checkonly = 0;
checkfirst = 0;
force = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'check',5)
        checkfirst =1;
    elseif strncmpi(varargin{j},'force',5)
        force =1;
    end
    j = j+1;
end

done = 0; 
return;
%doesn't quite work yet.  If this is called, see no XY plot so can't change
%it. Or can see misleading XY plot based on the wrong spikes file;
if DATA.state.online 
    if DATA.state.applylastcluster %for now ? use something else?
        done = 0;
        return;
    end
    for e = DATA.currentexpt(:)'
        p = DATA.probe;
        if isfield(DATA,'allcodes') && size(DATA.allcodes,1) >= e && size(DATA.allcodes,2) >= p
            if (checkonly || checkfirst) && ~isempty(DATA.allcodes{e,p}) && sum(DATA.allcodes{e,p})
                done(e) = 1;
            end
            if checkfirst 
                if length(done) >= e && done(e)
                    force = 1;
                else
                    force = 0;
                end
            end
                
            if length(DATA.allcodes{e,p}) == size(DATA.AllData.Spikes.codes,2) || force
                if force
                    DATA.AllData.Spikes.codes = [];
                end
                DATA.AllData.Spikes.codes(:,2) = DATA.allcodes{e,p};
                DATA.AllData.Spikes.times = DATA.alltimes{e,p};
                done(e) = 1;
            end
        end
    end
    done = sum(done);
end
if checkonly
    DATA = sum(done);
end