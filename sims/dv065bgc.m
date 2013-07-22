
    NC = 3200;
    x = [1:NC] - NC/2;
    WLD = 6.0;

    hw.toLens_mm   = 522;
    hw.screen_mm = [400 300]; 
    hw.screen_px = [1600 1200]; 
    
    PXD = WLD * tan(pi/180) * hw.toLens_mm * hw.screen_px(1) / hw.screen_mm(1);
    HDSP = round(PXD / 8.0);   % 90 deg phase disparity of fundamental +-45
    PXD = HDSP * 8.0;

    % FUNDAMENTAL
    S1l = 4/pi*sin(2*pi*(x+0.5+HDSP)/PXD); 
    S1r = 4/pi*sin(2*pi*(x+0.5-HDSP)/PXD); 
    % SW
    SW90.l = sign(S1l);
    SW90.r = sign(S1r);
    % MF
    MF90.l = SW90.l - S1l;
    MF90.r = SW90.r - S1r;

    PXD = WLD * tan(pi/180) * hw.toLens_mm * hw.screen_px(1) / hw.screen_mm(1);
    HDSP = round(PXD / 12.0);   % 60 deg phase disparity of fundamental +-30
    PXD = HDSP * 12.0;
    disp(PXD);
    % FUNDAMENTAL
    S1l = 4/pi*sin(2*pi*(x+0.5+HDSP)/PXD); 
    S1r = 4/pi*sin(2*pi*(x+0.5-HDSP)/PXD); 
    % 5F
    S5l = 4/(pi*5)*sin(2*pi*(x+0.5+HDSP)/(PXD / 5.0)); 
    S5r = 4/(pi*5)*sin(2*pi*(x+0.5-HDSP)/(PXD / 5.0)); 
    % 7F
    S7l = 4/(pi*7)*sin(2*pi*(x+0.5+HDSP)/(PXD / 7.0)); 
    S7r = 4/(pi*7)*sin(2*pi*(x+0.5-HDSP)/(PXD / 7.0)); 
    % MF
    MF60.l = sign(S1l) - S1l;
    MF60.r = sign(S1r) - S1r;
    % MF-5F
    MF5F.l = MF60.l - S5l;
    MF5F.r = MF60.r - S5r;
    % MF-5F-7F
    MF7F.l = MF5F.l - S7l;
    MF7F.r = MF5F.r - S7r;
    % MF phase scrambled
    for i=1:1,
        MFps{i}.l = 0;
        MFps{i}.r = 0;
        N = 3;
        while (PXD / N) > 2,
            rp = rand(1) * 2 * pi;
            MFps{i}.l = MFps{i}.l + 4/(pi*N)*sin(2*pi*(x+0.5+HDSP)/(PXD / N)+rp); 
            MFps{i}.r = MFps{i}.r + 4/(pi*N)*sin(2*pi*(x+0.5-HDSP)/(PXD / N)+rp); 
            N = N + 2;
        end;
%         figure;
         plot(MFps{i}.l, 'r');
         hold on;
         plot(MFps{i}.r((HDSP*2+1):end), 'b');
    end;
    

