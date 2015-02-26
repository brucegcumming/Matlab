function str =ShowString(DATA, Expt)
str = '';
if isfield(Expt.Stimvals,'rb') & Expt.Stimvals.rb ~= 0
    str = sprintf('rb%.0f',Expt.Stimvals.rb);
end
if sum(strcmp(Expt.Stimvals.et,{'Pp' 'Op'}))
    DATA.show.or = 1; %DATA not retured, so this is just local
end
if sum(strcmp(Expt.Stimvals.et,{'or'})) && Expt.Stimvals.st ==21
    DATA.show.sf = 1; %DATA not retured, so this is just local
    if sum(strcmp(Expt.Stimvals.e2,{'ob'}))
        DATA.show.or = 1; %DATA not retured, so this is just local
    end
end
fn = fields(DATA.show);

for k = 1:length(fn)
    if ~strcmp(fn{k},'times')
        if strcmp(fn{k},'Fr') & x > 1
            x = GetEval(Expt,fn{k});
            str = [str sprintf(' %s=%.1f',fn{k},x)];
        elseif DATA.show.(fn{k})
            x = GetEval(Expt,fn{k});
            if sum(strcmp(fn{k},{'jx' 'pi'}))
                str = [str sprintf(' %s=%.3f',fn{k},x)];
            elseif sum(strcmp(fn{k},{'or' 'sl' 'se' }))
                str = [str sprintf(' %s=%.0f',fn{k},x)];
            elseif strcmp(fn{k},'ntrials')
                str = [str sprintf(' %s=%.0f',fn{k},length(Expt.Trials))];
            elseif ischar(x)
                str = [str sprintf(' %s=%s',fn{k},x)];
            else
                str = [str sprintf(' %s=%.2f',fn{k},x)];
            end
        end
    end
end


