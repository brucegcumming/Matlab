function DATA = BuildGUI(DATA)

scrsz = get(0,'Screensize');
DATA.gui.scrsz = scrsz;
if isfield(DATA.gui,'FontSize')
    FontSize = DATA.gui.FontSize;
else
    FontSize = 0;
end

cw= scrsz(3)/140;
ch= scrsz(4)/60;
if scrsz(3) > 2000
    cw = 10;
elseif scrsz(3) > 1900
    cw = 10;
end
if scrsz(4) > 1200
    ch = 10;
elseif scrsz(4) == 1200
    ch = 8;
elseif scrsz(4) > 1100
    ch = 9;
else
    ch= 9;
end
if DATA.state.bigwindow(1) > 0
    cw = cw * DATA.state.bigwindow(1);
end
if DATA.state.bigwindow(2) > 0
    ch = ch * DATA.state.bigwindow(2);
end
DATA.plot.cw = cw;
DATA.plot.ch = ch;
DATA.user = GetUserName();
DATA.host = gethostname();
%fprintf('User is %s\n',DATA.user);
wsiz = [cw*40,ch*44];
SPACE = 0.01;
nch = 0.05;
ncw= 1./45; %45 chars wide

if length(DATA.layout.top) ~= 4
    DATA.layout.top = [100 scrsz(4)-(wsiz(2)+ch*4) wsiz(1) wsiz(2)];
end
cntrl_box = figure('Position', DATA.layout.top,...
    'NumberTitle', 'off', 'Tag',DATA.tag.top,'Name',DATA.tag.top,'ResizeFcn',@cmb.GuiResize);
DATA.toplevel = cntrl_box;

if ~isempty(DATA.configfile)
    prefconfig = [DATA.gui.prefsdir '/' DATA.configfile '.config'];
    if ~exist(DATA.configfile) && exist(prefconfig,'file')
        DATA.configfile = prefconfig;
    end
    DATA = cmb.DoConfig(DATA, DATA.configfile, 'load');
end

if DATA.options.uselastlayout && strcmp(DATA.layoutfile,DATA.defaultlayout)
    DATA = cmb.DoLayout(DATA, 'loadlast');
else
    DATA = cmb.DoLayout(DATA, DATA.layoutfile,'load');
end
if FontSize > 0
    DATA.gui.FontSize = FontSize;
elseif ~isfield(DATA.gui, 'FontSize')
    DATA.gui.FontSize = 12;
end

bp = [0.02 0.02 0.98 0.25];
%            set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);
set(cntrl_box,'DefaultUIControlFontSize',DATA.gui.FontSize);
DATA.elst = uicontrol(gcf, 'Style','listbox',...
    'Callback', ['cmb.combine(''setexpt'',''Tag'',''' DATA.tag.top ''');'],'Tag','subexptlist',...
    'Units','norm','Position',bp,'Max',3,'Min',1);
bp = [0.02 bp(2)+bp(4)+0.01 ncw*9 nch];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['cmb.combine(''combine'',''Tag'',''' DATA.tag.top ''');'],...
    'String', 'Combine','Units','norm', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = ncw * 7;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['cmb.combine(''Save'',''Tag'',''' DATA.tag.top ''');'],...
    'String', 'Save', 'Units','norm','Position', bp);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*7 nch];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['cmb.combine(''SaveLFP'',''Tag'',''' DATA.tag.top ''');'],...
    'String', 'SaveLFP', 'Units', 'norm', 'Position', bp);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*14 nch];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@cmb.SaveWithEM},...
    'String', 'Save with EM', 'Units', 'norm', 'Position', bp);


bp = [SPACE bp(2)+bp(4)+SPACE 0.99 nch];
DATA.saveitem = uicontrol(gcf,'Style','Edit','String','Save As','Units', 'norm','Position',bp);

bp = [0.01 bp(2)+bp(4)+SPACE 0.98 nch*6];
DATA.clst = uicontrol(gcf, 'Style','listbox',...
    'Callback', ['cmb.combine(''listexps'',''Tag'',''' DATA.tag.top ''');'],'Tag','explist',...
    'Units', 'norm','Position',bp,'Max',3,'Min',1);

bp = [0.01 bp(2)+bp(4)+SPACE ncw*9 nch];
uicontrol(gcf,'Style','Text','String','Data File','Units', 'norm','Position',bp);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*30 nch];
uicontrol(gcf,'Style','Edit','String',DATA.datafilename,'Units', 'norm','Position',bp,'Callback',['cmb.combine(''newfile''''Tag'',''' DATA.tag.top ''');'],...
    'Tag','FileName');

bp = [0.01 bp(2)+bp(4)+SPACE ncw*9 nch];
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'ShowN', 'Tag', 'ShowN', 'Units', 'norm','Position', bp, 'value',DATA.plot.showN);

bp(1) = bp(1) + bp(3) + SPACE;
bp(3) = ncw * 8;
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Spikes', 'Tag', 'QuickSpks', 'Units', 'norm','Position', bp,'value',DATA.state.showspikes,'Callback',@cmb.Update);
bp(1) = bp(1) + bp(3) + SPACE;
bp(3) = ncw * 9;
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'SpkXY', 'Tag', 'SpkXY', 'Units', 'norm','Position', bp,'value',DATA.state.showspkxy,...
    'Callback',@cmb.Update);

bp(1) = bp(1) + bp(3) + SPACE;
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Clean', 'Tag', 'ResetClusters', 'Units', 'norm','Position', bp,'value',DATA.state.resetclusters,...
    'Callback',@cmb.Update);

bp(1) = bp(1) + bp(3) + SPACE;
bp(3) = ncw*6;
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Auto', 'Tag', 'AutoPlot','Units', 'norm', 'Position', bp,'value',DATA.state.autoplot,...
    'Callback',@cmb.Update);

bp = [0.01 bp(2)+bp(4)+SPACE ncw*6 nch];
uicontrol(gcf,'Style', 'pop',...
    'String', 'Mean|Seq:time|Seq:Trials|Seq:Id|Seq:Edit|Psych Only|Collapse Y|Cp Trials|CP Hist|sdf|ACloop|ACresp|sdf for 0|Sdf collapse X|Sdf collapse Y|Counts|Subspace|Pcolor', ...
    'Tag', 'PlotSeq','Units', 'norm', 'Position', bp,'value',DATA.state.plotseq+1,...
    'Callback',@cmb.Update);

bp = [bp(1)+bp(3)+SPACE bp(2) ncw*9 nch];
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Recount', 'Tag', 'Recount', 'Units', 'norm','Position', bp,'value',DATA.state.recount,...
    'Callback',@cmb.Update);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*7 nch];
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Psych', 'Tag', 'PlotPsych', 'Units', 'norm','Position', bp,'value',DATA.state.plotpsych,...
    'Callback',@cmb.Update);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*11 nch];
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Combined', 'Tag', 'PlotCombined', 'Units', 'norm','Position', bp,'value',DATA.state.plotcombined,...
    'Callback',@cmb.Update);



bp = [0.01 bp(2)+bp(4)+SPACE ncw*3 nch];
uicontrol(gcf,'Style', 'checkbox',...
    'String', '0', 'Tag', 'UseCluster0', 'Units', 'norm','Position', bp,'Callback',{@cmb.SetClusters, DATA.tag.top});

bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
    'String', '1', 'Tag', 'UseCluster1', 'Units', 'norm','Position', bp,'Callback',{@cmb.SetClusters, DATA.tag.top});

bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
    'String', '2', 'Tag', 'UseCluster2', 'Units', 'norm','Position', bp,'Callback',{@cmb.SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
    'String', '3', 'Tag', 'UseCluster3', 'Units', 'norm','Position', bp,'Callback',{@cmb.SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
    'String', '4', 'Tag', 'UseCluster4', 'Units', 'norm','Position', bp,'Callback',{@cmb.SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3)+ncw/2;
bp(3) = ncw*4;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @cmb.FitButton, ...
    'String', 'Fit', 'Tag','FitButton','Units', 'norm','Position', bp);
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'cmb.combine(''setexpplot'')', ...
    'String', 'Plot', 'Tag','PlotButton','Units', 'norm','Position', bp);
bp(1) = bp(1)+bp(3);
bp(3) = ncw*6;
uicontrol(gcf,'Style', 'checkbox',...
    'String', 'LFP', 'Tag', 'UseLFP', 'Units', 'norm','Position', bp,'Callback',@cmb.Update);

bp(1) = bp(1)+bp(3);
h = uicontrol(gcf,'Style', 'pop','String',DATA.probenames,'Units', 'norm','Position', bp,...
    'Tag','ProbeId','Callback',@cmb.SetProbeHit,'value',1,'interruptible','off');
if length(DATA.probelist) > 1
    bp(1) = bp(1)+bp(3);
    uicontrol(gcf,'style','pushbutton','string','+','Position',cp, 'Callback',{@cmb.SetProbeHit,'next'});
end

bp = [0.01 bp(2)+bp(4)+SPACE 0.98 nch];
uicontrol(gcf,'Style','Edit','String','Comments','Units', 'norm','Position',bp,'callback',@cmb.AddComment);

DATA.gui.toprow = bp(2);

hm = uimenu(gcf,'Label','&Mark');
uimenu(hm,'Label','&Mark','Callback','cmb.combine(''Mark'')');
uimenu(hm,'Label','&relist','Callback','cmb.combine(''relist'')');
uimenu(hm,'Label','&Options','Callback','cmb.combine(''options'')');
uimenu(hm,'Label','&ShowVals','Callback','cmb.combine(''showvals'')');
uimenu(hm,'Label','&To Front','Callback','cmb.combine(''winfront'')');
uimenu(hm,'Label','&Update RF','Callback','cmb.combine(''rfupdate'')');
uimenu(hm,'Label','&Next CEell','Callback','cmb.combine(''nextcell'')');
uimenu(hm,'Label','&Update Lists','Callback','cmb.combine(''checklists'')');
uimenu(hm,'Label','&Comments','Callback','cmb.combine(''Comments'')');
%  uimenu(hm,'Label','Save Comments','Callback','cmb.combine(''SaveComments'')');
sm = uimenu(hm,'Label','&Settings');
uimenu(sm,'Label','Load &Config','Callback',{@cmb.DoConfig, 'select'});
uimenu(sm,'Label','&Save Config','Callback',{@cmb.DoConfig, 'save'});
uimenu(sm,'Label','&Save Default Config','Callback',{@cmb.DoConfig, 'savedefault'});
uimenu(sm,'Label','Load &Layout','Callback',{@cmb.DoLayout, 'choose'});
uimenu(sm,'Label','Save La&yout','Callback',{@cmb.DoLayout, 'save'});
uimenu(sm,'Label','Choose &Font','Callback',{@cmb.DoLayout, 'setfont'});
uimenu(sm,'Label','Save &DefaultLayout','Callback',{@cmb.DoLayout, 'savedefault'});
uimenu(hm,'Label','&Close, kill kept figures','Callback',{@cmb.CloseCombine, 'closeall'});
uimenu(hm,'Label','&Close','Callback',@cmb.CloseCombine);
cm = uimenu(gcf,'Label','&Cluster');
uimenu(cm,'Label','&AutoCut','Callback','cmb.combine(''autocut'')');
uimenu(cm,'Label','&AutoCluster','Callback','cmb.combine(''autocluster'')');
uimenu(cm,'Label','&AutoCluster All','Callback','cmb.combine(''allautocluster'')');
uimenu(cm,'Label','&TrackCluster','Callback','cmb.combine(''trackcluster'')');
uimenu(cm,'Label','&DDF','Callback','cmb.combine(''plotddf'')');
uimenu(cm,'Label','&FixMains','Callback','cmb.combine(''fixmains'')');
uimenu(cm,'Label','&ClearAll (one probe)','Callback','cmb.combine(''clearallclusters'')');
uimenu(cm,'Label','&ClearOnline','Callback','cmb.combine(''clearonlineclusters'')');
uimenu(cm,'Label','&Reload All','Callback',{@cmb.ReloadClusters, 'all'});
uimenu(cm,'Label','&Reload One Probe','Callback',{@cmb.ReloadClusters});
uimenu(cm,'Label','&Reload and Classify','Callback',{@cmb.ReloadClusters,'relcassify'});
uimenu(cm,'Label','&Make Templates','Callback',@cmb.MakeTemplates);
uimenu(cm,'Label','&Add Templates','Callback',@cmb.AddTemplate);
uimenu(cm,'Label','&Plot Templates','Callback',@cmb.PlotTemplates);
if DATA.state.online
    cm = uimenu(gcf,'Label','&Relist','Callback','cmb.combine(''relist'')');
end
cm = uimenu(gcf,'Label','&Options','Callback','cmb.combine(''options'')');
cm = uimenu(gcf,'Label','&Quick');
uimenu(cm,'Label','Online - Tuning only','Callback',{@cmb.OptionMenu,'QuickConfig'},'tag','onlinenospikes');
uimenu(cm,'Label','Online - with Spikes','Callback',{@cmb.OptionMenu,'QuickConfig'},'tag','onlinespikes');
uimenu(cm,'Label','Offline - Tuning only','Callback',{@cmb.OptionMenu,'QuickConfig'},'tag','offlinenospikes');
uimenu(cm,'Label','Offline - with Spikes','Callback',{@cmb.OptionMenu,'QuickConfig'},'tag','offlinespikes');

set(gcf,'Menubar','none');
DATA.toplevel = cntrl_box;
DATA.figpos(1).tag = DATA.tag.top;
DATA.figpos(1).pos = get(DATA.toplevel,'Position');
set(get(cntrl_box,'children'),'interruptible','off');
set(DATA.toplevel,'UserData',DATA);

setappdata(DATA.toplevel,'FigPos',DATA.figpos);

