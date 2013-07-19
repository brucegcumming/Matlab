%Expt = FillTrials(Expt, code, ...)
%Expt = FillTrials(Expt, code, 'mode') uses the modal value, not the mean
%FillTrials(Expt, 'saccade', [siz timeoffset timerange]) adds a
%vector 'PostSacc' to each Trial where 1 indicates stimuli that fell
%between timeoffset to timeoffset+timerange after a saccade. 
%
%
%
function Expt = FillTrials(Expt,  varargin)

j = 1;
while j <= length(varargin)
  if(strcmpi(varargin{j},'mode'))
    mode = 1;
  elseif(strncmpi(varargin{j},'consec',4))
      Expt = FillConsecScore(Expt,code);
  elseif(strncmpi(varargin{j},'saccade',4))
      j = j+1;
      Expt = FillSaccades(Expt, varargin{j});
  elseif strcmp(varargin{j},'BOPRC') %special case. Calculate dO/ce from
    Expt = FillBOPRC(Expt);
  else
    Expt = FillField(Expt, varargin{j});
  end
  j = j+1;
end

   
function Expt = FillConsecScore(Expt,type)

consec = 0;
if ~isfield(Expt.Trials,type)
    if ~isempty(Expt.Stimvals.e2)
        if isnumeric(Expt.Stimvals.e2)
            type = 'dx';
        else
            type = Expt.Stimvals.e2;
        end
    elseif isfield(Expt.Trials,'dx')
        type = 'dx';        
    end
end
if ~isfield(Expt.Trials,type)
    [Expt.Trials.consec] = deal(0);
    return;
end
scores = [Expt.Trials.RespDir] .* sign([Expt.Trials.(type)]);
Expt.Trials(1).consec = 0;
for j = 1:length(Expt.Trials)-1
    if scores(j) < 0
        consec = consec+1;
    elseif scores(j) > 0
        consec = 0;
    end
    Expt.Trials(j+1).consec = consec;
end


function Expt = FillField(Expt, code)

if strcmp(code,'sO') %fill this even if it exists already
%        if isfield(Expt.Trials,'yo') && isfield(Expt.Trials,'xo') && isfield(Expt.Trials,'or') && isfield(Expt.Trials,'wi') && isfield(Expt.Trials,'hi')
        if isfield(Expt.Trials,'yo') && isfield(Expt.Trials,'xo') && isfield(Expt.Trials,'or')
            rf = Expt.Stimvals.rf;
            [a,b] = Counts([Expt.Trials.or]);
            [c,d] = sort(a);
            ors = b(d);
            sins = sin(ors.*pi/180);
            coss = cos(ors.*pi/180);
            sina = sins(1);
            cosa = coss(1);
            if length(ors) == 4
                orthors = ors([3 4 1 2]);
            end
            Expt.Stimvals.xycos = min(abs([sina cosa]));
            for j = 1:length(Expt.Trials)
%                ar = Expt.Trials(j).wi./Expt.Trials(j).hi;
                ar = Expt.Trials(j).ar;
                oid = find(ors == Expt.Trials(j).or);
                sina = sins(oid);
                cosa = coss(oid);
                Expt.Trials(j).stO = [(rf(1) - Expt.Trials(j).xo) .* sina +  (Expt.Trials(j).yo-rf(2)) .* cosa] ;
                Expt.Trials(j).stP = [(rf(1) - Expt.Trials(j).xo) .* cosa +  (Expt.Trials(j).yo-rf(2)) .* sina] ;
                if (ar > 1 && abs(sins(oid)) > 0.9) || (ar < 1 && abs(sins(oid)) < 0.1)
                    Expt.Trials(j).sO = Expt.Trials(j).stO;
                    Expt.Trials(j).stimxy = Expt.Trials(j).xo;
                    Expt.Trials(j).xydir = 1; 
                elseif (ar > 1 && abs(sins(oid)) < 0.1) || (ar < 1 && abs(sins(oid)) > 0.9 )
                    Expt.Trials(j).sO = Expt.Trials(j).stP;
                    Expt.Trials(j).stimxy = Expt.Trials(j).yo;
                    Expt.Trials(j).xydir = 2;
                elseif (ar > 1 && sins(oid) < 0) || (ar < 1 && sins(oid) > 0 )
                    Expt.Trials(j).sO = Expt.Trials(j).stO;
                    Expt.Trials(j).stimxy = (Expt.Trials(j).xo+Expt.Trials(j).yo)./sqrt(2);
                    Expt.Trials(j).xydir = 1;
                elseif (ar > 1 && sins(oid) > 0) || (ar < 1 && sins(oid) < 0 )
                    Expt.Trials(j).sO = Expt.Trials(j).stO;
                    Expt.Trials(j).stimxy = Expt.Trials(j).yo;
                    Expt.Trials(j).stimxy = (Expt.Trials(j).xo-Expt.Trials(j).yo)./sqrt(2);
                    Expt.Trials(j).xydir = 2;
                else
                    Expt.Trials(j).stimxy = Expt.Trials(j).yo+rf(1)-rf(2) ;
                    Expt.Trials(j).sO = Expt.Trials(j).stP;
                    Expt.Trials(j).xydir = 5;
                end
                if ar > 1
                    Expt.Trials(j).thinor = Expt.Trials(j).or;
                else
                    Expt.Trials(j).thinor = orthors(oid);
                end
            end
        elseif isfield(Expt.Trials,'sO') %
            if ~isfield(Expt.Trials,'xo')
                for j = 1:length(Expt.Trials)
                    cosa = cos(Expt.Trials(j).or *pi/180);
                    sina = sin(Expt.Trials(j).or *pi/180);
                    Expt.Trials(j).xo = Expt.Trials(j).sO * cosa;
                    Expt.Trials(j).yo = Expt.Trials(j).sO * sina;
                    Expt.Trials(j).stimxy = Expt.Trials(j).sO * cosa;
                end
            end
        else
            for j = 1:length(Expt.Trials)
                Expt.Trials(j).sO = Expt.Trials(j).xo;
            end
        end
elseif ~isfield(Expt.Trials(1),code) 
    if strmatch(code,'oR')
        idx = find([Expt.Trials.Trial]);  
        for j = 1:length(idx)
            Expt.Trials(idx(j)).oR = Expt.Trials(idx(j)).or - Expt.Trials(idx(j)).od;
        end
    elseif strmatch(code,'rd') %relative disp
        idx = find([Expt.Trials.Trial]);  
        for j = 1:length(idx)
            Expt.Trials(idx(j)).rd = Expt.Trials(idx(j)).dx - Expt.Trials(idx(j)).bd;
        end
    elseif strmatch(code,'dur')
        idx = find([Expt.Trials.Trial]);  
        for j = 1:length(idx)
            Expt.Trials(idx(j)).dur = Expt.Trials(idx(j)).End(end) - Expt.Trials(idx(j)).Start(1);
        end
    elseif isfield(Expt.Stimvals,code) & ~isempty(Expt.Stimvals.(code))
        idx = find([Expt.Trials.Trial]);  
        eval(['[Expt.Trials(idx).' code '] = deal(Expt.Stimvals.' code ');']);
    elseif strcmp(code,'dO')
        if isfield(Expt.Trials,'dy')
        else
            for j = 1:length(Expt.Trials)
                Expt.Trials(j).dO = Expt.Trials(j).dx;
            end
        end
    elseif strmatch(code,'rndphase')
        idx = find([Expt.Trials.Trial]);
        for j = 1:length(idx)
            if isempty(strfind(Expt.Trials(idx(j)).OptionCode,'+rp'))
                Expt.Trials(idx(j)).rndphase = 0;
            else
                Expt.Trials(idx(j)).rndphase = 1;
            end2
            end
        
        end
    else %can happen with online files
        idx = find([Expt.Trials.Trial]);  
        [Expt.Trials.(code)] = deal(NaN);
    end
end


if isfield(Expt.Trials,code)
if length(cat(1,Expt.Trials.(code))) < length(Expt.Trials)
    val = eval(['mode([Expt.Trials.' code '])']);
    for j = 1:length(Expt.Trials)
        if eval(['isempty(Expt.Trials(j).' code ')']) & ~isempty([Expt.Trials(j).Start])
                eval(['Expt.Trials(j).' code ' = val;']);
        end
    end
end
end

function Expt = FillBOPRC(Expt)

for j = 1:length(Expt.Trials)
    Expt.Trials(j).dO = Expt.Trials(j).Ol-Expt.Trials(j).Or;
    id = find(Expt.Trials(j).Ol < -1000);
    Expt.Trials(j).dO(id) = Expt.Trials(j).Ol(id);
    ce = Expt.Trials(j).lph == Expt.Trials(j).rph;
    Expt.Trials(j).ce = (ce *2)-1;
end


function       Expt = FillSaccades(Expt, params);

minsz = params(1);
latency = params(2);
dur = params(3);

for j = 1:length(Expt.Trials)
    starts = Expt.Trials(j).Start;
    Expt.Trials(j).PostSacc = zeros(size(starts));
    Expt.Trials(j).Trigger = [];
    for k = 1:length(Expt.Trials(j).Saccades)
        S = Expt.Trials(j).Saccades(k);
        if S.size > minsz
        t(1) = Expt.Trials(j).TrialStart+S.start+latency;
        t(2) = t(1)+dur;
        si = find(starts > t(1) & starts <= t(2));
        Expt.Trials(j).PostSacc(si) = 1;
        Expt.Trials(j).Trigger(si) = S.start;
        end
    end
end

