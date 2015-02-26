function opdx = hv2op(hv, Ro)
%opdx = hv2op(hv, Ro)
%Given [hdisp vdisp] and Ro calculate dO and dP

if isstruct(hv) && isfield(hv,'Stimvals') %given an expt
    Expt = hv;
    for j = 1:length(Expt.Trials)
        if isfield(Expt.Trials,'dy')
            %don't overwrite subspace vectors
            for j = 1:length(Expt.Trials)
                op = hv2op([Expt.Trials(j).dx Expt.Trials(j).dy],Expt.Stimvals.Ro);
                if ~isfield(Expt.Trials,'dP') | isempty(Expt.Trials(j).dP) | isnan(Expt.Trials(j).dP)
                    Expt.Trials(j).dP = op(2);
                end
                if ~isfield(Expt.Trials,'dO') |  isempty(Expt.Trials(j).dO) | isnan(Expt.Trials(j).dO)
                    Expt.Trials(j).dO = op(1);
                end
            end
        end
    end
end
    ca = cos(Ro * pi/180);
    sa = sin(Ro * pi/180);
% need to recheck sign conventions here in replay...    

    opdx(1) = hv(1) .* sa - hv(2) .* ca;
    opdx(2) = hv(2) .* sa + hv(1) .* ca;
