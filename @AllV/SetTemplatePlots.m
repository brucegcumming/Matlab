function DATA = SetTemplatePlots(DATA, varargin)
%AllV sets the template spaces to plot, depending on probes, cluster
%version



    DATA.tmplots(1,:) = [1 3]; %p2(trigger) vs p1
    DATA.tmplots(2,:) = [1 4]; %p2(trigger) vs p3
    if DATA.cluster.version > 1.11
        if isfield(DATA,'TemplateScores') && size(DATA.TemplateScores,2) > 18
            DATA.tmplots(3,:) = [1 19]; %p2(trigger) vs p4 if present, else sum
        else
            DATA.tmplots(3,:) = [1 10]; %p2(trigger) vs p4 if present, else sum
        end
    else
        DATA.tmplots(3,:) = [1 10]; %p2(trigger) vs p4 if present, else sum
    end
    DATA.tmplots(4,:) = [2 10];
    DATA.tmplots(5,:) = [1 8];
    DATA.tmplots(6,:) = [2 11];
    DATA.tmplots(7,:) = [2 12];
    DATA.tmplots(8,:) = [1 5];
    DATA.tmplots(9,:) = [2 13];
    DATA.tmplots(10,:) = [2 14];
    DATA.tmplots(11,:) = [2 15];
    DATA.tmplots(12,:) = [1 6];
    DATA.tmplots(13,:) = [6 11];
    DATA.tmplots(14,:) = [2 16];
    DATA.tmplots(15,:) = [2 17];
    DATA.tmplots(16,:) = [2 18];

% tmplots 17:24 is for Std Template Plots    
    if isfield(DATA,'chspk')
        np = length(DATA.chspk);
    else
        np = 3;
    end
    nt = DATA.ntemplates*2;
    DATA.tmplots(17,:) = [nt+1 1]; %p2(trigger) vs p1 template 1
    DATA.tmplots(18,:) = [nt+1 2*nt+1]; %p2(trigger)T1 vs p3T1
    DATA.tmplots(19,:) = [nt+1+nt/2 1+nt/2]; %p2T2dt p1T1dt
    DATA.tmplots(20,:) = [nt+1+nt/2 2*nt+1+nt/2]; %p2T2dt p3T1dt
    DATA.tmplots(21,:) = [nt+1 nt+2]; %p2 1 vs2
    DATA.tmplots(22,:) = [nt+1 nt+3]; %p2 1 vs 3
    if nt > 3
        DATA.tmplots(23,:) = [nt+1 nt+4]; %p2 1 vs 4
        DATA.tmplots(24,:) = [nt+1+nt/2 nt+4+nt/2];
    else
        DATA.tmplots(23,:) = [nt+1+nt/2 2+nt/2];
        DATA.tmplots(24,:) = [3 4];
    end
    DATA.tmplots(25,:) = [nt+2 nt+3]; %p2 2 vs 3
    DATA.tmplots(26,:) = [nt+2 nt+4]; %p2 2 vs 4
    DATA.tmplots(27,:) = [nt+3 nt+4]; %p2 3 vs 4
    DATA.tmplots(28,:) = [nt+2+nt/2 nt+3+nt/2]; %p2 2 vs 3
    DATA.tmplots(29,:) = [nt+2+nt/2 nt+4+nt/2]; %p2 2 vs 4
    DATA.tmplots(30,:) = [nt+3+nt/2 nt+4+nt/2]; %p2 3 vs 4
    DATA.tmplots(31,:) = [nt+2 3]; %p2 2 vs p1 3
    DATA.tmplots(32,:) = [nt+2 nt*2+3]; %p2 2 vs p3 3
    DATA.tmplots(33,:) = [nt+2 2]; %p2 T2 vs P1T2
    DATA.tmplots(34,:) = [nt+2 nt*2+2]; %p2 T2 vs p3T2
    DATA.tmplots(35,:) = [nt+3 3]; %p2 T3 vs p3T3
    DATA.tmplots(36,:) = [nt+4 4]; %p2 T4 vs T4
    DATA.tmplots(37,:) = [nt+4 nt*2+4]; %p2 T4 vs p3T4
    DATA.tmplots(38,:) = [nt+2+nt/2 2+nt/2]; %p2 T2dt3 vs p1t2dt
    DATA.tmplots(39,:) = [nt+2+nt/2 (nt*2)+2+nt/2]; %p2 t2dt vs p3 t2dt
    DATA.tmplots(40,:) = [nt+2 nt*2+3]; %p2 2 vs p3 3

    if isfield(DATA,'TemplateScores')
        id = find(DATA.tmplots) > size(DATA.TemplateScores,2);
        DATA.tmpplots(id) = 1;
    end
    if DATA.nprobes == 1
        DATA.tmplots(4,:) = [1 10];
        DATA.tmplots(6,:) = [1 11];
        DATA.tmplots(7,:) = [1 12];
    end
