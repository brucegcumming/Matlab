function DATA = SetExpt(DATA, varargin)
%set fields in Expt neeeded in AllVPcs

Expt = DATA.Expt;
if isempty(Expt)
    return;
end
if ~isfield(Expt.Header,'expname')
    Expt.Header.expname = Expt2Name(Expt);
end

Expt.Header.trialdur = sum([Expt.Trials.dur]);
Expt = FillTrials(Expt,'ed');
Expt = FillTrials(Expt,'st');
[a,b] = fileparts(Expt.Header.Name);
Expt.Header.title = [b '.' Expt.Header.expname];
Expt.Header.preperiod = DATA.preperiod * 10000;
Expt.Header.postperiod = DATA.postperiod * 10000;
if ~isfield(Expt.Header,'ReadMethod')
    Expt.Header.ReadMethod = -1;
end
fprintf('Expt%d: %s %d trials\n',DATA.exptno,Expt.Header.title,length(Expt.Trials));
Expt.exptno = DATA.exptno;
Expt = FillTrials(Expt,'ed');
Expt = FillTrials(Expt,'st');

DATA.Expt = Expt;

