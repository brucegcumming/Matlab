function SetGUI(DATA)    SetCheck('SpikeXY',DATA.plot.spkxy,DATA.toplevel);    SetCheck('PlotTrigHist',DATA.plot.trighist,DATA.toplevel);    SetCheck('PlotHistogram',DATA.plot.hist,DATA.toplevel);    SetCheck('PlotXYseq',DATA.plot.xyseq,DATA.toplevel);    SetCheck('PlotSpks',DATA.plot.quickspks,DATA.toplevel);    SetCheck('LabelEd',DATA.show.ed,DATA.toplevel);    SetCheck('LabelExptName',DATA.show.exptname,DATA.toplevel);    SetCheck('LabelExptno',DATA.show.exptno,DATA.toplevel);    SetCheck('AllVPcs',DATA.show.allvpcs,DATA.toplevel);    SetCheck('SpikeMean',DATA.plot.spkmean,DATA.toplevel);    SetCheck('RefitGM',DATA.plot.refitgm,DATA.toplevel);    SetCheck('ShowGM',DATA.plot.gmfit,DATA.toplevel);        F = FindFig(DATA.tag.celllist);    f = fields(DATA.markcell);    onoff = {'off' 'on'};    for j = 1:length(f)        hm = findobj(F,'Tag',['mark' f{j}],'type','uimenu');        on = DATA.markcell.(f{j}) > 0;        set(hm,'checked',onoff{1+on});    end    f = {'watchallcellspks' 'watchallcellmean' 'watchallcellxy' 'watchexptallspks'};    for j = 1:length(f)        hm = findobj('Tag', f{j},'type','uimenu');        if ~isempty(hm)         on = DATA.show.(f{j}) > 0;        set(hm,'checked',onoff{1+on});        end    end        SetMenuCheck(DATA.toplevel, 'allvpcsmode', DATA.allvpcsmode,'exclusive');        it = findobj('Tag','DataShowMenu');    for k = 1:length(it)        c = get(it(k),'children');        for j = 1:length(c);            type = get(c(j),'tag');            if isfield(DATA.plot,type)                set(c(j),'checked',onoff{1+DATA.plot.(type)});            end        end    end