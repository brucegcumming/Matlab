function DATA = VergPysch(DATA, varargin)
%standalone version of Psych plotter from verg.

    
    
    
    if isempty(DATA.Expts) || isempty(DATA.Trials)
        return;
    end
    if strmatch(DATA.psych.blockmode,'None')
        return;
    end
    DATA = CheckExpts(DATA);

    Expt = [];
    e = length(DATA.Expts);
    allid = [];
    if strcmp(DATA.psych.blockmode,'All') || sum(DATA.plotexpts)
        if strmatch(DATA.psych.blockmode,'All')
            expts = 1:e-1;
        else
            expts = find(DATA.plotexpts);
        end
        %can only combine expts if have same et,e2 types.  Check each
        %against last in list
        if sum(strcmp(DATA.psych.blockmode,{'Current' 'All'}))
            Expt = DATA.Expts{e};
        else
            Expt = DATA.Expts{expts(end)};
        end
        
        for j = expts
            if isfield(DATA.Expts{j},'last')
                last = DATA.Expts{j}.last;
            else
                last = length(DATA.Trials);
            end                
            if strcmp(DATA.Expts{j}.Stimvals.et,Expt.Stimvals.et) && ...
               strcmp(DATA.Expts{j}.Stimvals.e2,Expt.Stimvals.e2)
                allid= [allid DATA.Expts{j}.first:last];
            end
        end
    else
        Expt = DATA.Expts{e};
    end
    %Always update trial list for current expt
    id = DATA.Expts{e}.first:length(DATA.Trials);
        DATA.Expts{e}.Trials = DATA.Trials(id);
    if Expt.first == DATA.Expts{e}.first  %using current
        if strmatch(DATA.psych.blockmode,{'All' 'Current' 'OnlyCurrent'})
         allid = [allid id];
        end
    end
    id = unique(allid);
    if length(id) < 2
        return;
    end
    Expt.Trials = DATA.Trials(id);
    Expt.Header.rc = 0;
    Expt.Header.expname  = 'Online';
    if isfield(Expt.Trials,'RespDir')
        Expt.Header.psych  = 1;
    else
        Expt.Header.psych  = 0;
    end
    
    f = {'m2' 'or' 'm3' 'n2'};
    for j = 1:length(f)
        if isfield(DATA.binoc{1},f{j})
            Expt.Stimvals.(f{j}) = DATA.binoc{1}.(f{j});
        end
    end    
    
    
    if DATA.psych.show
        DATA = SetFigure('VergPsych', DATA);
        hold off;
        eargs = {};
        if DATA.psych.trialresult
            plot([Expt.Trials.Start],[Expt.Trials.good],'o');
            datetick('x','hh:mm');
            axis('tight');
        elseif DATA.psych.showblocks
            ExptPsych(DATA.Expts(expts),'labelblock',eargs{:});
        else
        np = sum(abs([Expt.Trials.RespDir]) ==1); %psych trials
        Expt = FillTrials(Expt,Expt.Stimvals.et);
        Expt = FillTrials(Expt,Expt.Stimvals.e2);
        if np > 1
            if strcmp(Expt.Stimvals.e2,'od') && isfield(Expt.Stimvals,'or')
                for j = 1:length(DATA.expvals{2})
                    do = DATA.expvals{2}(j);
                    lo = Expt.Stimvals.or + do/2;
                    ro = Expt.Stimvals.or - do/2;
                    legendlabels{j} = sprintf('R%.0dL%.0f',ro,lo);
                end
                eargs = {eargs{:} 'legendlabels' legendlabels};
            end
            try
                [a,b] = ExptPsych(Expt,'nmin',1,'mintrials',2,'shown',eargs{:});
            catch ME
                fprintf('Error In ExptPsych %s\n',ME.message);
            end
            id = find(strcmp(Expt.Stimvals.et,{DATA.comcodes.code}));
            if length(id) == 1
                set(get(gca,'xlabel'),'string',DATA.comcodes(id).label);
            end
            if DATA.psych.crosshairs
                plot([0 0], get(gca,'ylim'),'k:')
                plot(get(gca,'xlim'),[0.5 0.5],'k:')
            end
        end
        end
    end
    DATA.Expt = Expt;
    
    function DATA = SetFigure(tag, DATA)

    [a,isnew] = GetFigure(tag);
    onoff = {'off' 'on'};
    if isnew
        DATA.figs.(tag) = a;
        if strcmp(tag,'VergPsych')
            hm = uimenu(a, 'Label','Expts','Tag','ExptMenu');
            PsychMenu(DATA);
            hm = uimenu(a, 'Label','Options','Tag','PsychOptions');
            sm = uimenu(hm,'Label','Crosshairs','callback', {@ChoosePsych, 'crosshairs'},...
                'checked',onoff{DATA.psych.crosshairs+1});
            sm = uimenu(hm,'Label','Just Show Trial outcomes','callback', {@ChoosePsych, 'trialresult'},...
                'checked',onoff{DATA.psych.trialresult+1});
            sm = uimenu(hm,'Label','Separate by Block','callback', {@ChoosePsych, 'showblocks'},...
                'checked',onoff{DATA.psych.showblocks+1});
            sm = uimenu(hm,'Label','Save Trial Data','callback', {@ChoosePsych, 'savetrials'});
            set(a,'UserData',DATA.toplevel);
            set(a,'DefaultUIControlFontSize',DATA.font.FontSize);
            set(a,'DefaultUIControlFontName',DATA.font.FontName);

        end
        set(DATA.toplevel,'UserData',DATA);
    end

function PsychMenu(DATA)    
    if ~isfield(DATA,'figs')
        return;
    end
    hm = findobj(DATA.figs.VergPsych,'tag','ExptMenu');
    c = get(hm,'children');
    delete(c);
    for j = 1:length(DATA.Expts)
        if isfield(DATA.Expts{j},'Trials')
        nt = length(DATA.Expts{j}.Trials);
        sm = uimenu(hm,'Label', sprintf('Expt%d %s %d',j,Expt2Name(DATA.Expts{j}),nt),'CallBack', {@ChoosePsych, sprintf('Expt%d',j)});
        if j < length(DATA.plotexpts) && DATA.plotexpts(j)
            set(sm,'Checked','on');
        end
        end
    end
    sm = uimenu(hm,'Label', 'Current','Callback',{@ChoosePsych, 'Current'});
    if strcmp(DATA.psych.blockmode,'Current')
        set(sm,'Checked','On');
    end
    uimenu(hm,'Label', 'Only Current','Callback',{@ChoosePsych, 'OnlyCurrent'},'tag','OnlyCurrent');
    uimenu(hm,'Label', 'All','Callback',{@ChoosePsych, 'All'},'tag','All');
    uimenu(hm,'Label', 'One','Callback',{@ChoosePsych, 'OneOnly'},'tag','OneOnly');
    uimenu(hm,'Label', 'None','Callback',{@ChoosePsych, 'None'},'tag','None');
    
function DATA = CheckExpts(DATA)

    for j = 1:length(DATA.Expts)
        if j < length(DATA.Expts) && ~isfield(DATA.Expts{j},'last')
            DATA.Expts{j}.last = DATA.Expts{j+1}.first -1;
        end
        DATA.Expts{j}.Header.exptno = j;
        DATA.Expts{j}.Header.rc = 0; %Cant' do rc online at the moment
    end
    
