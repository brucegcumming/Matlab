function Plot_PrefSF_LvR_CB(fig,j);
DEFINITIONS
handles = guidata(fig);
results = handles.results;
% See if this plot is already up for this object:
if handles.plotup(j)==1
    % if it is up, then clicking again deletes it
    set(handles.plot(j),'deletefcn','')
    delete(handles.plot(j))
    handles.plotup(j)=0;
    delete(handles.txt(j))
else
    % if the plot is not already up, cliucking brings it up
    Lbestfittype = results(j).SFTuning.fit.BestFitType(L);
    Lbestconstr = results(j).SFTuning.fit.BestConstraint(L);
    Rbestfittype = results(j).SFTuning.fit.BestFitType(R);
    Rbestconstr = results(j).SFTuning.fit.BestConstraint(R);
    SFL = results(j).SFTuning.fit.params(L,Lbestfittype,Lbestconstr).freqpref;
    SFR = results(j).SFTuning.fit.params(R,Rbestfittype,Rbestconstr).freqpref;
            
    disp(['Cell ' results(j).filebase ' : SFL = ' num2str(SFL) '; SFR = ' num2str(SFR) ]);
    handles.txt(j)=text(SFL,SFR,strrep(results(j).filebase,'.0',''),'fontsize',8,'verticalal','bot');
    
    handles.plot(j)=figure('pos',[739   551   424   357]);
    Plot_SFTC(results(j).SFTuning);
    title(results(j).filebase);
    set(handles.plot(j),'name',results(j).filebase)
    
    % record that it is up
    handles.plotup(j)=1;
    
    % Make sure that if the user closes the figure themselves, we keep track of it
    set(handles.plot(j),'deletefcn',['Plot_PrefSF_LvR_CB(' num2str(fig) ',' num2str(j) ')'])
end
guidata(fig,handles);