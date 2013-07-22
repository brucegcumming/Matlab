function Expt = FixExpt(Expt,type)

needpattern = 0;
if strcmp(type,'ed') & isfield(Expt.Trials,'ed')
    
    id = find(Expt.Trials.ed ~= 0);
    if length(id) > 1 & id(1) > 10 & diff(Expt.Trials.ed(id(1)-1:id(1))) > 0.1
        Expt.Trials.ed(1:id(1)) = Expt.Trials.ed(id(1));
    end
end

if strcmp(type,'dp')
    if isfield(Expt.Trials,'a2')
        needpattern = 1;
    end
    if needpattern
        for j = 1:length(Expt.Trials)
            or = Expt.Trials(j).or;
            a2 = or-Expt.Trials(j).a2;
            da = (or-a2) .* pi/180;
            angle = (or+a2) .*pi/180;
            dmag = Expt.Trials(j).dp .* cos(da);
            Expt.Trials(j).pdx = dmag .* cos(angle);
            Expt.Trials(j).pdy = dmag .* cos(angle);
        end
    end
end