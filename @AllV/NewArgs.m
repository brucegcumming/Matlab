function args = NewArgs(DATA,type)
% determine args required for a new call to AllVpcs
% when switching probe/expt
args = {};
if strcmp(type,'expt')
    if DATA.fullvswitchmode.applylast
        args = {args{:} 'applylast'};
    end   
    if DATA.fullvswitchmode.refcut
        args = {args{:} 'refcut'};
    end   
    if DATA.fullvswitchmode.summary
        args = {args{:} 'summary'};
    end   
end
