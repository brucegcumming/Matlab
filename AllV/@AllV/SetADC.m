function SetADC(a,b,fcn)    DATA = GetDataFromFig(a);    DATA.adcmousept.down = 0;    DATA.adcmousept.h = -1;            if fcn == 10 %callback when button released         DATA.vsmps(DATA.setadcpos) = DATA.adcmousept.x;        DATA.vpts = AllV.SetVsamples(DATA,DATA.probe, DATA.npallv, DATA.nvpts);        DATA.dvpts = DATA.vpts;        set(gcf,'ButtonDownFcn',{@AllV.ShowADCPos, 0});        set(gca,'ButtonDownFcn',{@AllV.ShowADCPos, 0});        set(gcf,'WindowButtonMotionFcn',{@AllV.ShowADCPos, 0});        set(gcf,'WindowButtonUpFcn',{@AllV.ShowADCPos, 0});        DATA.plottype = 2;        SetData(DATA);        AllV.CheckOneXY(DATA);        AllV.ReplotPCs(DATA,[]);            else        DATA.setadcpos = fcn;        set(gcf,'ButtonDownFcn',{@AllV.ShowADCPos, 1});        set(gca,'ButtonDownFcn',{@AllV.ShowADCPos, 1});        set(gcf,'WindowButtonMotionFcn',{@AllV.ShowADCPos, 2});        set(gcf,'WindowButtonUpFcn',{@AllV.ShowADCPos, 3});    end    set(DATA.toplevel,'UserData',DATA);    