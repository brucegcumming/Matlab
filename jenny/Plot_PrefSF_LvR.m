% Gets up various plots I am interested in
function PrefSF_LvR(results)
% Plot peak SF in L eye against peak SF in R eye
cnt=0;
DEFINITIONS

figure('pos',[10   137   755   657],'name','Frequency scatterplot','numbertitle','off')

%Mark on identity
plot([-1 1]*100,[-1 1]*100,'k')
% this has to be under the plot, otherwise it interferes with the user interaction

for j=1:length(results)
    disp([num2str(j) ' : '  results(j).filebase])
    if isstruct(results(j).SFTuning) % this results was available for this cell
        try % becase if not all data were available, will get an error here

            Lbestfittype = results(j).SFTuning.fit.BestFitType(L);
            Lbestconstr = results(j).SFTuning.fit.BestConstraint(L);
            Rbestfittype = results(j).SFTuning.fit.BestFitType(R);
            Rbestconstr = results(j).SFTuning.fit.BestConstraint(R);
            SFL = results(j).SFTuning.fit.params(L,Lbestfittype,Lbestconstr).freqpref;
            SFR = results(j).SFTuning.fit.params(R,Rbestfittype,Rbestconstr).freqpref;
            SFLsd = results(j).SFTuning.fit.params(L,Lbestfittype,Lbestconstr).SD;
            SFRsd = results(j).SFTuning.fit.params(R,Rbestfittype,Rbestconstr).SD;
            % This may be empty if the results was not recorded monocularlu
            if ~isempty(SFR) & ~isempty(SFL)
                cnt=cnt+1;
                SF_L(cnt) = SFL;
                SF_R(cnt) = SFR;
                SF_Lsd(cnt) = SFLsd;
                SF_Rsd(cnt) = SFRsd;
                % If fit was linear, the SD is in freq:
                hold on
                if Lbestfittype == LIN
                    % do the error bars as +/- SD of fitted Gaussian
                    plot(SFL+[1 -1]*SFLsd,SFR*[1 1],'c');
                elseif Lbestfittype == LOG
                    % Then it's a bit more complicated as the fit was done in log space
                    % Plot error bars as +/- 1 SD of fitted Gaussian in log space
                    SFLmin = SFL / exp(SFLsd);%=exp(log(SFL) - SFLsd);
                    SFLmax = SFL * exp(SFLsd);
                    plot([SFLmin SFLmax],SFR*[1 1],'c');
                end
                if Rbestfittype == LIN
                    plot(SFL*[1 1],SFR+[1 -1]*SFRsd,'c');
                elseif Rbestfittype == LOG
                    % Then it's a bit more complicated as the fit was done in log space
                    % Plot error bars as +/- 1 SD of fitted Gaussian in log space
                    SFRmin = SFR / exp(SFRsd);
                    SFRmax = SFR * exp(SFRsd);
                    plot(SFL*[1 1],[SFRmin SFRmax],'c');
                end

                SFplotcallback{cnt} = ['Plot_PrefSF_LvR_CB(gcf,' num2str(j) ')'];
                
            end
        catch
        end % try/catch
    end % if isstruct
end
                
% Store the fact that no figures have presently been called up for this figure:
handles = guihandles(gcf);
for j=1:length(results)
    handles.plotup(j)=0;
end
handles.results = results;
guidata(gcf,handles);
% Make sure that, when you close this figure, you close all the little windows that popped up from it
set(gcf,'deletefcn',['CloseChildWindows(' num2str(gcf) ')'])

xlabel('Preferred freq in L eye (cpd)')
ylabel('Preferred freq in R eye (cpd)')
title(['Scatterplot of preferred SF in the two eyes; ' num2str(cnt) ' cells'])
set(gca,'xlim',[-0.5 max(SF_L)*1.1],'ylim',[-0.5 max(SF_R)*1.1]);
%mark on identity
plot(get(gca,'xlim'),get(gca,'xlim'),'k')
plot(get(gca,'ylim'),get(gca,'ylim'),'k')

% Finally plot the blobs marking the data-points. These have to go last otherwise the buttondown clalback doesn't work
for j=1:cnt
    blob(j) = plot(SF_L(j),SF_R(j),'bo','markerfacecolor','b','buttondownfcn',SFplotcallback{j},'markersize',5);
end
