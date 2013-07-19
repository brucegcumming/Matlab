function res = runnsine(im, varargin)
%
%runnsine(im, ....)
%Analyze/Plot Nsine RC data 
%
%if im is the image strucure recording all the components for each see
%then this builds a set of sdfs for each condition.
%
% runnsine(Expname,'build') makes all the files required for the analysis
% based on then name in Expt.Header.
% runnsine(Exptname,'save') makes all the files and saves them.
%
%if im is a set of SDFs build on a previous call, this plots them up
%different ways.
%
% Plot types:
%  'monoc', eye
%       plots monocular responses to SF. Eye is a number determining what data
%       plotted - 1 = Left, 2 =  Right, 0 = Sum of L+R, 3 = 2 subplots, L and R
%                 4 = 3 subplot L, R and Binoc (collapsed over phase)
% 'corrs', corrs
%           where corrs is a matrix of binocular xcorrelation functions (1
%           for each image. Plots the spike trigggered cross-correlation.
%  'cmpcorr', corrs    
%           plots the spike-triggered cross-correlation and the predicted
%           value based on summing the dp tuning functions (fit with sines)
%           for each frequency component.
%  'cmpf', corrs    
%           plots the FT of the spike-triggered cross-correlation and the binocular SF tuning;
%
% 'sinamps'  plots amplitudes of DP modulation for each frequency wrt time.
%
%
%If arg 2 is an Expt file, this is used, otherwise it is loaded based on
%the file named in im.name.  If you have already loaded and Expt, its
%quicker if you pass it as an argument.
%


%
% Todo - make plot for comparing binocular mean and binocular amps
% for these comparison plots doing slices if named, else pcolors, add
% an autoslice option. - If 2 maps, be sure and cover max var for both.
% compare SF tuning, L, R, L+R and binoc at differnt slices.
%
% Make a plot that checks all dtfiles plots with im.dtfile plot, showing n;
DPPLOT = 1;
DPSFIMS = 2;
STC = 3;
DISP = 4;
PLOTALL = 5;
DPSFSLICES = 6;
SFDPIM = 7;
SFPLOT = 8;
XCORR = 9;
XCORRFT = 10;
DPVARS = 11;
MONOCSF = 12;
DTPREDICT = 13;
DTPREDICTA = 14;
SINAMPS = 15;
CMPCORR = 16;
CMPF = 17;
CMPA = 18;
CMPSF = 19;  %plots SF tuning for revcor and SDFs 
ACPLOT = 20;
SUMMARY = 21;
STA = 22;
nsdfs = [];
Expt = [];
res = [];
plottype = DPPLOT;
meanmode = 0;
slices = [20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95];
xyslices = [];
xarg = {};
buildall = 0;
printfigs = 0;
savesummary = 0;
labela = [];
sfslices = [];
recalculate = 0;
eyelabels = {'L+R','Left','Right'};
impath = '/home/mat1/bgc/Images/stimlists/';
if ispc
    if exist('Z:/files/STC','dir')
        stcprefix = 'Z:/files/STC';
    else
        stcprefix = 'C:/bgc/files/STC';
    end
    impath = 'X:/bgc/Images/stimlists/';
else
    stcprefix = '/bgc/files/STC';
    impath = '/home/mat1/bgc/Images/stimlists/';
end



if isstruct(im) & isfield(im,'rcres') & isfield(im,'im') % a summary file
    rcres = im.rcres;
    xc = im.xcdat;
    corrs = xc.corrs;
    nsdfs = im.nsdfs;
    sumdat = rmfield(im,{'rcres','xcdat', 'nsdfs'});
    im = im.im;
    sumdat = rmfield(sumdat,'im');
end

if isstruct(im) & isfield(im,'rcres') & isfield(im,'ftc') % a summary file witout raw data
    rcres = im.rcres;
    sumdat = rmfield(im,{'rcres'});
    clear im;
    im.name = sumdat.name;
    if ~isfield(sumdat,'pixsz')
        sumdat.pixsz = 0.1;
        sumdat.npix = 32;
    else
        im.xvals = sumdat.pixsz * ([1:sumdat.npix] - sumdat.npix/2);
%        im.lpixels = im.xvals;
    end
  
    xc.disps = rcres.corrx;
end

if isfield(im,'allr')
    nsdfs = im;
end


colors = {};
j = 1;
while j < nargin
    if isstruct(varargin{j})
        if isfield(varargin{j},'Trials')
            Expt = varargin{j};
        elseif isfield(varargin{j},'allr')
            nsdfs = varargin{j};
        elseif isfield(varargin{j},'rcid')
            rcres = varargin{j};
            im.delays = rcres.delays;
            im.rcid = rcres.rcid;
            if isfield(rcres,'sfs') & ~isfield(im,'sfvals')
                im.sfvals = rcres.sfs;
            end
            if isfield(rcres,'dpvals') & ~isfield(im,'dpvals')
                im.dpvals = rcres.dpvals;
            end
        elseif isfield(varargin{j},'corrs') & isfield(varargin{j},'disps')
            xc = varargin{j};
            corrs = xc.corrs;
        elseif isfield(varargin{j},'rcres') & isfield(varargin{j},'im') % a summary file
            im = varargin{j}.im;
            rcres = varargin{j}.rcres;
            xc = varargin{j}.xcdat;
            corrs = xc.corrs;
            nsdfs = varargin{j}.nsdfs;
        end
    elseif strncmpi(varargin{j},'ac',2)
        plottype = ACPLOT;
    elseif strncmpi(varargin{j},'disp',3)
        plottype = DISP;
    elseif strncmpi(varargin{j},'all',3)
        plottype = PLOTALL;
    elseif strncmpi(varargin{j},'cmpf',4)
        plottype = CMPF;
        if j < nargin-1 & isstruct(varargin{j+1})
            j = j+1;
            corrs  = varargin{j};
        end
    elseif strncmpi(varargin{j},'cmpa',4)
        j = j+1;
        if ~exist('xc','var')
            corrs  = varargin{j};
        end
        plottype = CMPA;
    elseif strncmpi(varargin{j},'cmpsf',4)
        plottype = CMPSF;
    elseif strncmpi(varargin{j},'cmpcorr',4)
        j = j+1;
        corrs  = varargin{j};
        plottype = CMPCORR;
    elseif strncmpi(varargin{j},'corft',5)
        j = j+1;
        corrs  = varargin{j};
        plottype = XCORRFT;
    elseif strncmpi(varargin{j},'corrs',3)
        plottype = XCORR;
        if j < nargin-1 & isstruct(varargin{j+1})
            j = j+1;
            corrs  = varargin{j};
        end
    elseif strncmpi(varargin{j},'delay',3)
        j = j+1;
        xyslices  = varargin{j};
    elseif strncmpi(varargin{j},'dtpredict',3)
        plottype = DTPREDICT;
    elseif strncmpi(varargin{j},'dtfitpredict',3)
        plottype = DTPREDICTA;        
    elseif strncmpi(varargin{j},'dpsf',3)
        plottype = DPSFIMS;
    elseif strncmpi(varargin{j},'dpv',3)
        plottype = DPVARS;
    elseif strncmpi(varargin{j},'eigs',3)
        j = j+1;
        eigs = varargin{j};
    elseif strncmpi(varargin{j},'label',3)
        j = j+1;
        labela = varargin{j};
    elseif strncmpi(varargin{j},'monoc',3)
        plottype = MONOCSF;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            eye = varargin{j};
        else
            eye = 0;
        end
    elseif strncmpi(varargin{j},'nomean',3)
        meanmode = 1;
    elseif strncmpi(varargin{j},'plot',3)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'print',3)
        printfigs = 1;
    elseif strncmpi(varargin{j},'sfdp',4)
        plottype = DPSFSLICES;
    elseif strncmpi(varargin{j},'sfslices',4)
        j = j+1;
        sfslices = varargin{j};
    elseif strncmpi(varargin{j},'sf',2)
        plottype = SFPLOT;
    elseif strncmpi(varargin{j},'sinamps',3)
        plottype = SINAMPS;
    elseif strncmpi(varargin{j},'rc',2)
        j = j+1;
        rcres = varargin{j};
        im.delays = rcres.delays;
        im.rcid = rcres.rcid;
        if isfield(rcres,'sfs') & ~isfield(im,'sfvals')
            im.sfvals = rcres.sfs;
        end
        if isfield(rcres,'dpvals') & ~isfield(im,'dpvals')
            im.dpvals = rcres.dpvals;
        end
    elseif strncmpi(varargin{j},'recalc',4)
        recalculate = 1;
    elseif strncmpi(varargin{j},'slices',3)
        j = j+1;
        slices = varargin{j};
    elseif strncmpi(varargin{j},'sta',3)
        plottype = STA;
    elseif strncmpi(varargin{j},'stc',3)
        plottype = STC;
    elseif strncmpi(varargin{j},'sumsave',4)
        plottype = SUMMARY;
        savesummary = 1;
    elseif strncmpi(varargin{j},'summary',3)
        plottype = SUMMARY;
    elseif strncmpi(varargin{j},'times',3)
        j = j+1;
        xarg = {xarg{:} 'times' varargin{j}};
    end
    j = j+1;
end

if isempty(colors)
    colors = mycolors;
end
madeimages = 0;
madexcorr = 0;
if isstr(im)
    if strncmpi(im,'build',2)
        buildall = 1;
    elseif strncmpi(im,'save',2)
        buildall = 2;
    elseif strfind(im,'files/STC')
        load(im);
    else  %%named a cell number. load files
        cellname = im;
        cellname = strrep(cellname,'.c1','');
        cellname = strrep(cellname,'.c2','');
        cid = regexp(im,'\.c[0-9]');
        if ~isempty(cid)
            clusternum = str2num(im(cid(1)+2));
        else
            clusternum = 1;
        end
        imfile = sprintf('%s/%snsines.mat',stcprefix,cellname);
        spkfile = sprintf('%s/%s.c%d.spks.mat',stcprefix,cellname,clusternum);
        xcfile = sprintf('%s/%scorrs.mat',stcprefix,cellname);
        sdffile = sprintf('%s/%s.c%d.sdfs.mat',stcprefix,cellname,clusternum);
% once nsres was saved in im. This is not a good idea when omore than one
% cluster is associated with an image sequence.
        if exist(imfile,'file') & ~recalculate & clusternum ==1 & ~exist(spkfile)
            clear nsres;
            load(imfile);
            if exist('nsres','var')
                nsres.name = im.name;
                nsres.dtfile = im.dtfile;
                nsres.dtfiles = im.dtfiles;
                nsres.pddpfile = im.pddpfile;
                save(spkfile,'nsres');
            end
        end
        if (exist(imfile,'file') & exist(spkfile,'file')) & ~recalculate
            clear nsres;
            load(imfile);
            load(spkfile);
            res.rcres = nsres;
            rcres = nsres;
            rcres.nspikes = length(nsres.rcid{1});
            if isfield(nsres,'name')
                im.name = nsres.name
            end
            im.delays = rcres.delays;
            im.rcid = rcres.rcid;
            res.im = im;
            clear nsres;
            if exist(sdffile,'file')
                load(sdffile);
                nsdfs = nsdf;
                clear nsdf;
                res.nsdfs = nsdfs;
            end
            if exist(xcfile,'file')
                load(xcfile);
                res.xcdat = xc;
            end
        else
            load(name2path(sprintf('%s.c%d.nsines.DP.mat',cellname,clusternum)));
            if exist(imfile,'file') & ~recalculate  %% images already made - just redo spikes
                madeimages = 1;
            end
            if nargin > 1 & isstr(varargin{1})
                if strncmpi(varargin{1},'build',2)
                    buildall = 1;
                elseif strncmpi(varargin{1},'save',2)
                    buildall = 2;
                end
            end
        end
    end
end

if buildall
    if isempty(Expt)
        fprintf('Need Expt structure:\nrunnsine(''build'', Expt..)\n');
        return;
    end
    cellname = splitpath(Expt.Header.Name,'cellpref');
    clusternum = str2num(Expt.Cluster(1));
    %read in text table describing stimulus frames
    if madeimages
        load(imfile);
    else
        All = dlmread(sprintf('%s%sXnsines.txt',impath,cellname));
        %convert this to a useable structure
        im = tab2im(All,Expt);
    end
    im.name = splitpath(Expt.Header.Name);
    nsres = NsineRC(Expt, im, slices * 10, {});
	im.dpvals = unique(im.phases);
    nsres.dpvals = im.dpvals;
    nsres.name = im.name;
    if exist(xcfile,'file') & ~recalculate
        load(xcfile); 
        res.xcdat = xc;
        madexcorr = 1;
    elseif isfield(im,'lsfs')   %% Nsines random components for STC
        fprintf('Building Images from Fourier Components\n');
        xc  = CalcNsineCorr(im,'image');  %Calculate images
        im.lpixels = xc.lpixels;
        im.rpixels = xc.rpixels;
        im.xvals = xc.xvals; %locations of resulting pixels, in degrees.
        fprintf('Calculating cross-correlation (normalized)\n');
        dmax = floor((size(xc.lpixels,1)-1)/2);
        xc.corrs  = CalcDispCorr(im,[1:size(xc.lpixels,1)],-dmax:dmax,'circular');
        xc.disps = [-dmax:dmax] * mean(diff(im.xvals));
        nsres.corrx = xc.disps;
    else
        xc = CalcNsineCorr(im);
    end
    fprintf('Buidling SDF...');
    nsdf = BuildSdfs(im, Expt, res, xarg);
    if isfield(nsdf,'buildtime')
        fprintf(' Took %.1f sec',nsdf.buildtime);
    end
    fprintf('\n');
    tic;
    fprintf('Buidling Corr SDF...');
    corrsdf = Nsinesdf(Expt, im, 'rc', nsres, 'corrs', xc.corrs, 'box', 104);
    corrsdf.corrx = xc.disps;
    fprintf(' Took %.1f sec',toc);

    nsres = FindDTFiles(nsres);
    if buildall == 2 %build and save
        if ~madeimages
            fprintf('Saving %s\n',imfile);
            save(imfile,'im');
        end
        fprintf('Saving %s\n',spkfile);
        save(spkfile,'nsres');
        nam = sprintf('%s/%s.c%d.sdfs.mat',stcprefix,cellname,clusternum);
        fprintf('Saving %s\n',sdffile);
        save(sdffile,'nsdf','corrsdf','nsres');
        if ~madexcorr
            fprintf('Saving %s\n',xcfile);
            save(xcfile,'xc');
        end
    end
    res.im = im;
    res.corrsdf = corrsdf;
    res.nsres = nsres;
    res.nsdf = nsdf;
    rcres = nsres;
    rcres.nspikes = length(nsres.rcid{1});
    im.delays = rcres.delays;
    im.rcid = rcres.rcid;
    nsdfs = nsdf;
    clear nsdf;
    clear nsres;
   
    if plottype ~= SUMMARY
        return;
    end
end

if exist('xc','var')
    corrx = xc.disps;
    rcres.corrx = xc.disps;
elseif exist('rcres','var')
    corrx = rcres.corrx;
end

if exist('rcres','var')
    if isfield(rcres,'delays');
        im.delays = rcres.delays;
        if isfield(rcres,'rcid')
            im.rcid = rcres.rcid;
        end
    end
    if isfield(rcres,'sfs') & ~isfield(im,'sfvals')
        im.sfvals = rcres.sfs;
    end
    if isfield(rcres,'dpvals') & ~isfield(im,'dpvals')
        im.dpvals = rcres.dpvals;
    end
end

if plottype == SUMMARY
    % need to add binocular amps to F-T plots, and vars to time-var plot

    if ~isfield(res.rcres,'dtfile')
        res.rcres = FindDTFiles(res.rcres);
    end

    GetFigure('NsineRCa');
    cb = [];
    if exist('sumdat','var') & ~recalculate
        res = sumdat;
        subplot(3,2,6);
        fillpmesh(im.delays/10, im.sfvals, abs(res.predxc.sincmps),'plot');
        title('DP Amp');
        cb = colorbar;
    else
        xc.corrs(find(isnan(xc.corrs))) = 0;
        ftc = runnsine(im, rcres, xc, 'cmpf');  %% FT of ST-Xcorr vs binoc means
        xcr = runnsine(im, rcres, xc, 'corr');  %% ST-xcorr
        pcr = runnsine(im, rcres, xc, nsdfs, 'dtfit'); %% predict xcorr from dpsf
        res.ftc = ftc;
        res.predxc.im = pcr.im;
        res.predxc.disps = pcr.corrx;
        res.predxc.sincmps = pcr.sincmps;
        res.xc = xcr.xc;
        res.xcvar = std(xcr.xc');
        if isfield(im,'lsfs')
            res.mcx = runnsine(im, rcres, xc, 'monoc');
        end
        res.nframes = length(im.seeds);
        res.rcres = rcres;
        if isfield(im,'xvals')
            res.pixsz = mean(diff(im.xvals));
        elseif isfield(rcres,'corrx')
            res.pixsz = mean(diff(rcres.corrx));
        end
        if exist('nsdfs','var')
            bcx = runnsine(im, rcres, xc, nsdfs, 'cmpa');
            subplot(3,2,6);
            fillpmesh(im.delays/10, im.sfvals, abs(bcx.sincmps),'plot');
            res.dpamp = mean(abs(bcx.sincmps));
            title('DP Amp');
            cb = colorbar;
        end
    end
    SetColorBar(cb);

    subplot(3,2,1);
    
    hold off;
    fillpmesh(im.delays/10, im.sfvals, res.ftc.cft','plot');
    title('FT of xcorr');

    subplot(3,2,2);
    hold off;
    fillpmesh(im.delays/10, im.sfvals, res.ftc.sfim','plot');
    title('Binoc SF');

    subplot(3,2,4);
    hold off;
    if range(xc.disps) > 1/min(im.sfvals);
        id = find(abs(xc.disps) < 0.5/min(im.sfvals));
        imagesc(im.delays/10, xc.disps(id), res.xc(:,id)');
    else
        imagesc(im.delays/10, xc.disps, res.xc');
    end
    title('Xcorr');
    ylabel('disparity');
    SetColorBar(colorbar);

    if isfield(res,'mcx')
        dominance = max(res.mcx.lsfim(:)-0.5)/max(res.mcx.rsfim(:)-0.5);
        if dominance > 1
            subplot(3,2,3);
            hold off;
            fillpmesh(im.delays/10, im.sfvals, res.mcx.lsfim','plot');
%            imagesc(im.delays/10, im.sfvals, res.mcx.lsfim');
            title(sprintf('L (dominant %.1f) SF',dominance));
            
            subplot(3,2,5);
            fillpmesh(im.delays/10, im.sfvals, res.mcx.rsfim','plot');
            title('R SF');
        else
            subplot(3,2,5);
            hold off;
            fillpmesh(im.delays/10, im.sfvals, res.mcx.lsfim','plot');
            title('L SF');
            
            subplot(3,2,3);
            fillpmesh(im.delays/10, im.sfvals, res.mcx.rsfim','plot');
            title(sprintf('R (dominant %.1f) SF',1/dominance));
        end
    else
        subplot(3,2,3);
        ClearPlot(gca);
        subplot(3,2,3);
        ClearPlot(gca);
        if ~isfield(im,'lpixels')
            subplot(3,2,3);
            id = find(abs(res.predxc.disps) < 0.5/min(im.sfvals));
            imagesc(im.delays,res.predxc.disps(id), res.predxc.im(:,id)');
        end
    end
    
    
    
    subplot(3,2,1);
    text(0.8,1.2,sprintf('%s %d spikes %d frames',im.name,rcres.nspikes,res.nframes),'units','norm')
    subplot(3,2,3);
    text(0-0.2,0.5,'SF (cpd)','units','norm','rotation',90);
    subplot(3,2,5);
    text(1,-0.2,'Delay (ms)','units','norm');
    SetSubplots(3,2,[3 5],'caxis');

    GetFigure('NsineRCb');
    hold off;
    if exist('sumdat','var') & ~recalculate
        if isfield(res,'bestdelay')
            best = res.bestdelay;
        else
            best = round(length(rcres.delays)/2);
        end
        stcs.delaytimes = rcres.delays./10;
    elseif isfield(im, 'lpixels') %% If can't make pixels, can't do STC
        stcs = PlotSTC(im, rcres, xc, 'cmat','label','NsineRCb','noplot');
        res.stcvar = stcs.stcvar;
        res.stavar = stcs.stavar;
        res.staim = stcs.staim;
        res.xcvar = std(xcr.xc');
        best = stcs.bestdelay(1);
        res.bestdelay = best;
        res.name = im.name;
        if isfield(res.rcres,'dtfile')
            res.dtfile = res.rcres.dtfile;
            res.dtfiles = res.rcres.dtfiles;
            res.pddpfile = res.rcres.pddpfile;
        elseif isfield(im,'dtfile')
            res.dtfile = im.dtfile;
            res.dtfiles = im.dtfiles;
            res.pddpfile = im.pddpfile;
        end
        res.npix = size(im.lpixels,1);
    end

    if isfield(rcres,'bestdelay')
        best = rcres.bestdelay;
    end
    subplot(3,2,1);
    hold off;
    if isfield(res,'stcvar'); %% can do STC
        plot(stcs.delaytimes,res.stcvar);
        
        stavars = res.stavar .* max(res.stcvar(:))/max(res.stavar(:));
        hold on;
        plot(stcs.delaytimes,stavars,'--');
        plot(stcs.delaytimes,res.xcvar .* max(res.stcvar(:))/max(res.xcvar),'r--');
        legend('STC','L','R','B','STA','Xcorr');
        ylabel('Var');
        xlabel('Delay');
        
        if ~exist('sumdat') | recalculate
            stcres = PlotSTC(im, rcres, xc, 'cmat','delayid', best,'noplot');
            res.stcim = stcres.stcim;
        end
        
        subplot(3,2,2);
        hold off;
        imagesc(res.stcim - diag(diag(res.stcim)));
        title(sprintf('Covariance at %.0f',stcs.delaytimes(res.bestdelay)));
        hold on;
        nv = res.npix * 2;
        plot([0 nv],[0.5 + nv/2 0.5 + nv/2],'w');
        plot([0.5 + nv/2 0.5 + nv/2],[0 nv],'w');
        plot(0.5+[0 nv/2],0.5+[nv/2 nv],'w--');
        colorbar;
        xlabel(sprintf('Pixel (%.3f deg each)',res.pixsz));
        
        subplot(3,2,3);
        hold off;
        imagesc(stcs.delaytimes,im.xvals*2,res.staim');
        title('STA');
        xlabel('ms');
        
        %
        %  N.B. This plots disp tuning in subplot(2,3,5) also
        %
        subplot(3,2,4);
        hold off;
        if ~exist('sumdat') | recalculate
            dxres = PlotSTC(im, rcres, xc, 'eigx',[1 2 nv nv-1],'delayid', best,'subd',[3 2 4; 3 2 5]);
        else
            dxres = PlotSTC(im, rcres, sumdat, 'eigx',[1 2 nv nv-1],'delayid', best,'subd',[3 2 4; 3 2 5]);
        end
        if  ~isfield(dxres,'dxscale') % no DT data
            dxres.dxscale = [0 1];
        end
        xlabel('disparity');
        
    
        subplot(3,2,5);  %% Just add prediction fomr sine components
        if isfield(res,'predxc') %% predicted from sin component DP tuning
            pdxs = res.predxc.im(best,:);
            pdxs = dxres.dxscale(1) + pdxs * dxres.dxscale(2)/range(pdxs);
            predh = plot(res.predxc.disps, pdxs,'m');
        end
        if isfield(dxres,'dxdata')
            legend([dxres.dxdata.legend.h(1) predh dxres.handles(1:3)],{'odx' 'sfDP' 'diag' 'ev' 'xc'});
            legend('boxoff');
        end
        axis tight;
    else
        % need to make STA at least for this conditions
        ClearPlot(subplot(3,2,2));
        ClearPlot(subplot(3,2,5));
        ClearPlot(subplot(3,2,4));
        if ~exist('sumdat') | recalculate
            if ~isfield(im,'xvals')
                step = 1/(min(im.sfvals) * (2 * length(im.sfvals)+1));
                sz = 0.5/min(im.sfvals);
                xm = sz - step/2;
                im.xvals = -xm:step:xm; %%1 period, at Nyquist limit, 0 at center
            end
            
            for j = 1:length(im.delays)
                res.staim(j,:) = CalcSTA(im, im.rcid{j},im.xvals);
            end
            res.npix = length(im.xvals);
            res.stavar = std(res.staim');
        end
        subplot(3,2,1);
        hold off;
        plot(im.delays/10,res.dpamp);
        hold on;
        rscale = std(res.dpamp)./std(std(res.xc'));
        plot(im.delays/10,std(res.xc') * rscale,'r');
        rscale = std(res.dpamp)./std(std(res.staim'));
        plot(im.delays/10,std(res.staim') * rscale,'g');
        [amax, best] = max(std(res.xc'));
        legend('DP amp','XC');
        subplot(3,2,5);
        if isfield(res.rcres,'dtfile');
            tres = PlotExpt(name2path(im.nsres.dtfile));
            tvar = max(std(tres.means));
            tmean = mean(tres.means(:));
            trange = range(tres.means(:,end)); %% just use corr if there are both.
            id = find(xc.disps >= min(tres.x(:)) & xc.disps <= max(tres.x(:)));
            rscale = trange/range(res.xc(best,id));
            plot(xc.disps(id),tmean + res.xc(best,id) * rscale,'g')
            pdxs = res.predxc.im(best,:);
            pdxs = tmean + pdxs * trange/range(pdxs);
            id = find(res.predxc.disps >= min(tres.x(:)) & res.predxc.disps <= max(tres.x(:)));
            predh = plot(res.predxc.disps(id), pdxs(id),'m');
            axis tight;
        end
        subplot(3,2,4);
        [peak, peakf] = max(res.ftc.cft(best,:)); 
        hold off; plot(im.sfvals,res.ftc.cft(best,:));
        xfreqs = [1:size(res.ftc.xc,2)]/range(rcres.corrx);
        hold on; 
        pwsp =abs(fft(res.ftc.xc(best,:)));
        fid = find(xfreqs <= max(im.sfvals));
        
        plot(xfreqs(fid),pwsp(fid) * peak/max(pwsp),'g');
        pfreqs = [1:size(res.predxc.im,2)]/range(res.predxc.disps);
        fid = find(pfreqs <= max(im.sfvals));
        pwsp = abs(fft(res.predxc.im(best,:)));
        plot(pfreqs(fid),pwsp(fid) * peak/max(pwsp),'r');
        legend('cft','xc','dpred');
        
        
    end
    
    subplot(3,2,6);
    if ~exist('sumdat') | recalculate
        ac = runnsine(im, rcres, xc, nsdfs,'ac','delay',best);
        res.ac = ac;
    end
    hold off;
    if isfield(res.rcres,'pddpfile')
        pddp = PlotExpt(res.rcres.pddpfile,'legendpos',7);
        hold on;
        [mdiff, cid] = min(abs(pddp.linevals));
        [mdiff, aid] = min(abs(pddp.linevals-pi));
        rscale = range(pddp.means(:,cid))./range(res.ac.drf(:,1));
        id = find(xc.disps >= min(pddp.x(:)) & xc.disps <= max(pddp.x(:)));
        plot(xc.disps(id),mean(pddp.means(:)) + res.ac.drf(id,:).*rscale,'linewidth',2);
        legend(pddp.handles([cid aid]),{'C','A'});
    else
        plot(xc.disps,res.ac.drf);
    end
    subplot(3,2,1);
    text(0.8,1.2,sprintf('%s %d spikes %d frames',im.name,rcres.nspikes,res.nframes),'units','norm')
    
    GetFigure('NsineRCc');
    subplot(3,2,1);
    if exist('dxres','var')
        plot(dxres.eigvals);
        nv = length(dxres.eigvals);
        lr = polyfit([1:nv-8],dxres.eigvals(5:end-4)',1);
        rande = lr(2) + [1:nv] *lr(1);
        eigdev = abs(dxres.eigvals - rande');
        [a, res.eigrank] = sort(eigdev,'descend');
        for j = [1 3 5 7]
            text(res.eigrank(1),dxres.eigvals(res.eigrank(j)),sprintf('%d',j));
        end
        res.eigdev = eigdev;
        res.eigvals = dxres.eigvals;
        res.eig = dxres.eig;
        PlotSTC(im, rcres, res, 'eigval','delayid', best,'subid',[3 2 1]);
    title('EigenValues');
    a = 1; b =2;
    xl = get(gca,'Xlim');
    PlotSTC(im, rcres, res, 'eiglrs',[a:b],'delayid', best,'subid',[3 2 2]);
    legend('L','R');
    title(sprintf('Eigenvectors %d and %d',a,b));
    a = 3; b =4;
    PlotSTC(im, rcres, res, 'eiglr',[a:b],'delayid', best,'subid',[3 2 3]);
    title(sprintf('Eigenvectors %d and %d',a,b));
    b = length(dxres.eigvals);
    a = b-1;
    PlotSTC(im, rcres, res, 'eiglr',[a:b],'delayid', best,'subid',[3 2 4]);
    title(sprintf('Eigenvectors %d and %d',a,b));
   
    end
    subplot(3,2,5);
    hold off;
    npix = size(res.staim,2)/2;
    lid = 1:npix;
    rid = lid+npix;
    [amax, abest] = max(res.stavar);
    sta = res.staim(abest,:);
    plot(im.xvals,sta(lid),'-');
    hold on;
    plot(im.xvals,sta(rid),'--');
    staxc = xcorr(sta(rid),sta(lid),floor(length(lid)/2));
    rscale = std(sta)./std(staxc);
    plot(im.xvals,staxc * rscale,'r');
    set(gca,'xlim',[min(im.xvals) max(im.xvals)]);
    legend('L','R','XC');
    title(sprintf('STA at %.1fms',im.delays(abest)/10));
    
    if isfield(res.rcres,'sffile')
        subplot(3,2,6);
        PlotExpt(name2path(res.sffile));
        legend('boxoff');
    end
    
    subplot(3,2,1);
    text(0.8,1.2,sprintf('%s Eigenvectors at %.1f ms',im.name,im.delays(best)/10),'units','norm')
    
    if printfigs
        GetFigure('NsineRCa');
        set(gcf,'PaperPositionMode','manual','PaperUnits','normalized','PaperPosition',[0 0 1 1]);
        print;
        GetFigure('NsineRCb');
        set(gcf,'PaperPositionMode','manual','PaperUnits','normalized','PaperPosition',[0 0 1 1]);
        print;
        GetFigure('NsineRCc');
        set(gcf,'PaperPositionMode','manual','PaperUnits','normalized','PaperPosition',[0 0 1 1]);
        print;
    end
    if savesummary
        sumdat = res;
        sumdat.name = im.name;
        sumdat.xvals = im.xvals;
        if isfield(sumdat,'xcdat') sumdat = rmfield(sumdat,'xcdat'); end
        if isfield(sumdat,'im') sumdat = rmfield(sumdat,'im'); end
        if isfield(sumdat,'nsdfs') sumdat = rmfield(sumdat,'nsdfs'); end
        if isfield(sumdat,'rcres') sumdat.rcres = rmfield(sumdat.rcres,'rcid'); end
        nm = sprintf('/bgc/bgc/anal/rc/stcdat/%s', strrep(im.name,'DP','DPsum'));
        save(nm,'sumdat');
    end

    return;
elseif plottype == PLOTALL
    mcx = runnsine(im, 'monoc', 'label','Monoc SF');
    %need to put data into returned struct here
    %
    if ~isfield(res.rcres,'dtfile')
        res.rcres = FindDTFiles(res.rcres);
    end
    res.lsfim = mcx.lsfim;
    res.rsfim = mcx.rsfim;
    mcx = runnsine(im, rcres, xc, 'cmpf','label', 'Compare FT/F');
    res.binocsf = mcx.sfim;
    res.ftim = mcx.cft;
    res.corrim = mcx.xc;
    GetFigure('ST xcorr');
    subplot(1,1,1);
    mcx = runnsine(im, rcres,  xc, 'corr','label','ST xcorr');
    stcs = PlotSTC(im, rcres, xc, 'cmat','label','STCS');
    res.stcvar = stcs.stcvar;
    res.stavar = stcs.stavar;
    res.staim = stcs.staim;
    best = rcres.delays(stcs.bestdelay(1))/10;
    res.bestdelay = best;
    stcres = PlotSTC(im, rcres, xc, 'cmat','delays', best, 'label', 'BestSTC');
    res.stcim = stcres.stcim;
    mcx = PlotSTC(im, rcres, xc, 'eiglr',[1:2],'delays', best, 'label', 'Eigenvectors');
    res.evals = mcx.eigvals;
    res.evecs = mcx.eig;
    GetFigure('DTuning');
    subplot(1,1,1);
    mcx = PlotSTC(im, rcres, xc, 'eigx',[1:2],'delays', best, 'label', 'Eigenvector xcorr');
    res.sfvals = rcres.sfs;

    return;
end
    
if plottype == STA
    for j = 1:length(im.rcid)
        res.staim(j,:) = CalcSTA(im, im.rcid{j}, im.xvals);
    end
end

if ismember(plottype, [ACPLOT])
    if isfield(rcres,'bestdelay')
        delay = rcres.bestdelay;
    elseif ~isempty(xyslices)
        delay = xyslices(1);
    else
        delay = 7;
    end
    cths = prctile(corrs,[0 20 40 60 80 100]);
    for j = 1:length(xc.disps)
        cvs = corrs(im.rcid{delay},:);
        for k = 1:size(cths,1)-1
            res.ac(j,k) = sum(cvs(:,j) > cths(k,j) & cvs(:,j) <= cths(k+1,j))/sum(corrs(:,j) > cths(k,j) & corrs(:,j) <= cths(k+1,j));
        end
    end
    plot(res.ac(:,5));
    hold on;
    plot(res.ac(:,4),'m');
    plot(res.ac(:,3),'g');
    plot(res.ac(:,2),'y');
    plot(res.ac(:,1),'r');
    res.drf = [res.ac(:,5) - res.ac(:,3) res.ac(:,1) - res.ac(:,3)];
    pim = length(im.rcid{delay})./size(corrs,1);
    plot([1 size(res.ac,1)],[pim pim],':k');
    legend('Corr','0.5', '0', '-0.5', 'AC');
    return;
end
if ismember(plottype, [XCORR XCORRFT CMPCORR CMPF CMPA]) %% all need xcorr data
    for j = 1:length(im.delays)
        res.xc(j,:) = mean(corrs(im.rcid{j},:));
        res.ftxc(j,:) = abs(fft(res.xc(j,:)));
        for k = 1:length(im.sfvals)
            res.cft(j,k) = famp(corrx,res.xc(j,:),im.sfvals(k));
        end
    end
    GetFigure(labela);
    if ismember(plottype,[CMPCORR CMPF CMPA]) & isempty(xyslices)
        subplot(1,2,1);
    end
    if ismember(plottype,[XCORRFT CMPF CMPA])
        hold off;
        if isempty(xyslices)
            [X,Y] = meshgrid( im.delays);
            [X,Y,Z] = fillpmesh(im.delays./10,im.sfvals, res.cft');
            pcolor(X,Y,Z);
            shading('flat');
        else
            for j = 1:length(xyslices)
                [a, id] = min(abs(xyslices(j) - im.delays));
                plot(im.sfvals,res.cft(id,:),'color',colors{j});
                xymeans(j) = mean(res.cft(id,:));
                xystd(j) = std(res.cft(id,:));
                hold on;
            end
        end
        figure(gcf);
        title('FT of xcorr');
    elseif plottype == XCORRFT
        hold off;
        fc = 2:20;
        freqs = im.sfvals(1) .* (fc-1);
        imagesc(freqs, im.delays, res.ftxc(:,2:20));
        figure(gcf);
    else
        hold off;
        idx = find(abs(rcres.corrx) < 1); 
        imagesc(rcres.corrx(idx), im.delays./10, res.xc(:,idx));
        figure(gcf);
    end
    if ~ismember(plottype,[CMPCORR CMPF CMPA]);
        return;
    end
end

if ~isfield(nsdfs,'allr') & ismember(plottype,[SINAMPS DPSFSLICES DPSFIMS])
    fprintf('Building sdfs...');
    nsdfs = BuildSdfs(im, Expt, res, xarg);
    res.nsdfs = nsdfs;
    fprintf('\n');
end

if isfield(nsdfs,'allr')  %% this is a results file, not the data file
% allr is a 3-D array: allr (d,f,t) is the response at each sdf is a 2-D image of SF vs time, for a particular dp
    if ~isempty(labela)
        GetFigure(labela);
    end
    for j = 1:length(slices)
        [a, islices(j)] = min(abs(slices(j) - nsdfs.time./10));
    end
    
    if plottype == DPSFIMS
        [nr, nc] = Nsubplots(size(nsdfs.allr,1)); %dps
        for j = 1:size(nsdfs.allr,1)
            subplot(nr,nc,j);
            imagesc(nsdfs.time,im.sfvals,squeeze(nsdfs.allr(j,:,:)));
            title(sprintf('dp %.0f',im.dpvals(j)));
        end
        SetSubplots(nr,nc,1:j,'caxis');
    elseif ismember(plottype,[SFPLOT CMPF CMPSF])
        res.im = squeeze(mean(nsdfs.allr,1));
        if ismember(plottype,[CMPF]) & isempty(xyslices)
            subplot(1,2,2);
            hold off;
        elseif ismember(plottype,[CMPSF]) & isempty(sfslices)
            subplot(1,2,1);
            hold off;
        else
            subplot(1,1,1);
            if plottype == CMPF
                hold on;
            end            
        end
        if isempty(xyslices) & isempty(sfslices)
            [X,Y,Z] = fillpmesh(nsdfs.time./10, im.sfvals, res.im);
            pcolor(X,Y,Z);
            shading('flat');
            colorbar;
        elseif isempty(sfslices)
            for j = 1:length(xyslices)
                [a, id] = min(abs(xyslices(j) - nsdfs.time));
                sfdat = res.im(find(~isnan(res.im(:,id))),id);
                sfscale = xystd(j)/std(sfdat);
                scaled = sfscale * (res.im(:,id) - mean(sfdat)) + xymeans(j);
                plot(nsdfs.sfvals,scaled,':','color',colors{j});
                hold on;
            end
        else
            hold off;
            for j = 1:length(sfslices);
                sfr(j,:) = mean(squeeze(nsdfs.allr(:,sfslices(j),:)));
                plot(nsdfs.time,sfr(j,:),'color',colors{j});
                hold on;
            end
        end

    elseif ismember(plottype,[DPSFSLICES DPVARS])
        [nr, nc] = Nsubplots(length(slices));
        for j = 1:length(slices)
            dpsfim{j} = squeeze(nsdfs.allr(:,:,islices(j)));
        end
        if plottype == DPSFSLICES
            for j = 1:length(slices)-1
                subplot(nr,nc,j);
                imagesc(im.dpvals,im.sfvals,dpsfim{j});
                title(sprintf('delay  %.0f',slices(j)));
            end
            subplot(nr,nc,j+1);
            imagesc(im.dpvals,im.sfvals,nsdfs.netspikes);
            title(sprintf('Net Spikes',slices(j)));
            
        SetSubplots(nr,nc,1:length(slices)-1,'caxis');
        elseif plottype == DPVARS
            xticks = exp(linspace(log(0.5), log(16),6))
            xticks = xticks(find(xticks >= min(nsdfs.sfvals) & xticks <= max(nsdfs.sfvals)))
            xsc = 'log';
            subplot(2,1,1);
            hold off;
            for j = 1:length(slices);
%                plot(nsdfs.dpvals,std(dpsfim{j}));
                Z(:,j) = std(dpsfim{j});
                plot(nsdfs.sfvals,std(dpsfim{j}),'color',colors{j});
                hold on;
            end
            set(gca,'xscale',xsc,'Xtick',xticks);
            subplot(2,1,2);
            hold off;
            [X,Y,Z] = fillpmesh(slices,nsdfs.sfvals,Z);
            pcolor(Y,X,Z);
            shading('flat');
            colorbar;
            set(gca,'xscale',xsc,'Xtick',xticks);
        end
    elseif plottype == SFDPIM %DP vs time, one SF each panel
        [nr, nc] = Nsubplots(length(nsdfs.sfvals));
        for dp = 1:length(nsdfs.sdfs)
            for j = 1:size(nsdfs.sdfs{dp},2) %% # of time samples
                for k = 1:size(nsdfs.sdfs{dp},1)
                    dpsfim{k}(dp,j) = nsdfs.sdfs{dp}(k,j);
                end
            end
        end
        for k = 1:length(nsdfs.sfvals)
            subplot(nr,nc,k);
            imagesc(nsdfs.time,nsdfs.dpvals,dpsfim{k});
            title(sprintf('sf %.2f',nsdfs.sfvals(k)));
        end
        SetSubplots(nr,nc,1:k,'caxis');
    elseif plottype == MONOCSF
        if eye == 3
            subplot(2,1,1);
            res.rm = squeeze(nsdfs.monoc(2,:,:));
            res.lm = squeeze(nsdfs.monoc(1,:,:));
            if meanmode
                bmean = repmat(mean(res.rm + res.lm),size(nsdfs.monoc,2),1);
                res.rm = res.rm-bmean;
                res.lm = res.lm-bmean;
            end
            imagesc(res.rm);
            subplot(2,1,2);
            imagesc(res.lm);
            SetSubplots(2,1,[1:2],'caxis');
        elseif eye == 4
            
            subplot(3,1,1);
            res.rm = squeeze(nsdfs.monoc(2,:,:));
            res.lm = squeeze(nsdfs.monoc(1,:,:));
            if meanmode
                bmean = repmat(mean(res.rm + res.lm),size(nsdfs.monoc,2),1);
                res.rm = res.rm-bmean;
                res.lm = res.lm-bmean;
            end
            imagesc(res.rm);
            subplot(3,1,2);
            imagesc(res.lm);
            subplot(3,1,3);
            
            res.im = squeeze(mean(nsdfs.allr,1));
            if isempty(xyslices)
                imagesc(nsdfs.time./10, im.sfvals, res.im);
            end
            SetSubplots(3,1,[1:3],'caxis');
            colorbar;
        else
            subplot(1,1,1);
            if eye == 0 %% sum
                monoc = squeeze(mean(nsdfs.monoc,1));
            else
                monoc = squeeze(nsdfs.monoc(eye,:,:));
            end
            if meanmode
                rmean = repmat(mean(monoc),size(nsdfs.monoc,2),1);
                res.im = monoc-rmean;
            else
                res.im = monoc;
            end
            imagesc(res.im);
            title(eyelabels{eye+1});
            colorbar;
        end
    elseif ismember(plottype,[DTPREDICT DTPREDICTA SINAMPS CMPCORR CMPA])
        if ismember(plottype,[CMPCORR CMPA])
            subplot(1,2,2);
            hold off;
        else
            subplot(1,1,1);
        end
        if ~isfield(res,'corrx')
            res.corrx = linspace(-1,1,100);
        end
        if isfield(im,'dpvals')
            dpvals = [im.dpvals; min(im.dpvals) + 2 * pi];
        elseif isfield(rcres,'dpvals');
            dpvals = [rcres.dpvals; min(rcres.dpvals) + 2 * pi];
        else
            dpvals = [0:pi/6:pi];
        end
                
        edp = length(dpvals);
        allr = nsdfs.allr(:,:,islices);
        for j = 1:length(slices)
            for k = 1:size(nsdfs.allr,2)
                allr(edp,k,j) = allr(1,k,j);
            end
        end
        for j = 1:length(im.sfvals)
            dpivals(j,:) =  mod(2 * pi * res.corrx * im.sfvals(j), 2 * pi);
        end
        for j = 1:length(slices)
            ivals = [];
            for k = 1:size(nsdfs.allr,2);
                idx  = find(~isnan(allr(:,k,j)));
                if plottype == DTPREDICT
                    ivals(k,:) = interp1(dpvals(idx),squeeze(allr(idx,k,j)),dpivals(k,:));
                    idx = find(isnan(ivals(k,:)));
                    if ~isempty(idx)
                        ivals(k,idx) = mean(ivals(k,find(~isnan(ivals(k,:)))));
                    end
                    ptype = 'from interpolation';
                else
                    [a,c] = famp(dpvals(idx),squeeze(allr(idx,k,j)),1/(2 * pi));
                    if(length(idx) < 4)
                        a = 0;
                        c = 0;
                    end
                    ivals(k,:) = a * cos(-dpivals(k,:) + angle(c));
                    res.sincmps(k,j) = c;
                    ptype = 'from sin dp amps';
                end
            end
            res.im(j,:) = mean(ivals);
        end
        res.slices = slices;
        if ismember(plottype,[SINAMPS CMPA])
            [X,Y,Z] = fillpmesh(slices,im.sfvals,abs(res.sincmps));
            pcolor(X,Y,Z);
            title('Amplitude of DP resp');
            shading('flat');
        else
            imagesc(res.corrx,slices,res.im);
            title(sprintf('Predicted Corr resp %s',ptype));
        end
    end
    figure(gcf)
    if ~ismember(plottype,[CMPSF])
        return;
    end
end

if isempty(Expt)
    load(name2path(im.name));
end

if plottype == DPPLOT
    nsdf = Nsinesdf(Expt, im,'dp',3,'sdfw',20);
elseif plottype == DPSFIMS
    [nr,nc] = Nsubplots(length(im.dpvals));
    if ~isempty(nsdfs)
        res = nsdfs;
    else
        res = BuildSdfs(im, Expt, res, xarg)
    end
    for j = 1:length(im.dpvals);
        subplot(nr,nc,j);
        imagesc(squeeze(res.allr(j,:,:)));
        caxis(res.crange);
        title(sprintf('dp=%.2f',res.allr(j,1,1)));
    end
    colorbar;
elseif ismember(plottype,[MONOCSF])
    for j = 1:length(im.delays)
        lpc(j,:) = mean(im.lsfs(im.rcid{j},:));
        rpc(j,:) = mean(im.rsfs(im.rcid{j},:));
    end
    subplot(2,1,1);
    [X,Y,Z] = fillpmesh(im.delays./10, im.sfvals,lpc');
    pcolor(X,Y,Z);
    shading('flat');
    title('Left');

    subplot(2,1,2);
    [X,Y,Z] = fillpmesh(im.delays./10, im.sfvals,rpc');
    pcolor(X,Y,Z);
    shading('flat');
    title('Right');
    res.lsfim = lpc;
    res.rsfim = rpc;

elseif ismember(plottype,[SFPLOT CMPF CMPSF])
    for j = 1:length(im.delays)
        spc(j,:) = mean(im.sfs(im.rcid{j},:));
    end
    res.sfim = spc;
    if ~isempty(sfslices)
        for j = 1:length(sfslices);
            sj = sfslices(j);
            resp = (spc(:,sj) - mean(spc(:,sj))) .* std(sfr(j,:))./std(spc(:,sj));
            plot(im.delays,resp+mean(sfr(j,:)),':','color',colors{j});
            hold on;
        end
    elseif ismember(plottype,[CMPF CMPSF])
        subplot(1,2,2);
        hold off;
    else
        subplot(1,1,1);
    end
    if isempty(sfslices)
        [X,Y,Z] = fillpmesh(im.delays./10, im.sfvals,spc');
        pcolor(X,Y,Z);
        shading('flat');
    end
end
res.sfvals = im.sfvals;
res.dpvals = im.dpvals;


function res = BuildSdfs(im, Expt, res, xarg)

crange = [NaN NaN];


msdf = Nsinesdf(Expt, im,'monocs','box',104,'vals',im.sfvals,'image',xarg{:});
if isfield(msdf,'lsdf')
    for k = 1:length(im.sfvals)
        for ij = 1:size(msdf.rsdf,2)
            res.monoc(1,k,ij) = msdf.lsdf(k,ij);
            res.monoc(2,k,ij) = msdf.rsdf(k,ij);
            res.mn = msdf.n;
        end
    end
end
res.buildtime = 0;
for j = 1:length(im.dpvals);
    nsdf = Nsinesdf(Expt, im,'sfdp',im.dpvals(j),'box',104,'vals',im.sfvals,'image',xarg{:});
    res.buildtime = res.buildtime + nsdf.buildtime;
    cx = caxis;
    crange(1) = min([cx(1) crange(1)]);
    crange(2) = max([cx(2) crange(2)]);
    sdfs{j} = nsdf.sdfs;
    res.time = nsdf.times;
    res.n(j,:) = nsdf.n;
end
for j = 1:length(im.dpvals)
    for k = 1:length(im.sfvals)
        for ij = 1:size(sdfs{1},2)
            res.allr(j,k,ij) = sdfs{j}(k,ij);
            res.netspikes(j,k) = sum(sdfs{j}(k,:) - nsdf.mean);
        end
    end
end
res.crange = crange;
res.name = im.name;

function im = FindDTFiles(im)

ns = name2path(im.name);
pddp = strrep(ns,'DP','PDDP');
if exist(pddp,'file')
   im.pddpfile = pddp;
end
n = 1;
a = strrep(ns,'nsines.DP','rls.OXAC');
if exist(a,'file') im.dtfiles{n} = a; n = n+1; end
a = strrep(ns,'nsines.DP','rds.OXAC');
if exist(a,'file') im.dtfiles{n} = a; n = n+1;  end
a = strrep(ns,'nsines.DP','rls.ODX');
if exist(a,'file') im.dtfiles{n} = a;  n = n+1; end
a = strrep(ns,'nsines.DP','rds.ODX');
if exist(a,'file') im.dtfiles{n} = a;  n = n+1; end
a = strrep(ns,'nsines.DP','rds.AC');
if exist(a,'file') im.dtfiles{n} = a;  n = n+1; end
a = strrep(ns,'nsines.DP','rds.DT');
if exist(a,'file') im.dtfiles{n} = a;  n = n+1; end

%
% later, get this to load files and check if one has more n, etc.
% NB this all assumes user has build matlab files;
if n > 1
    im.dtfile = splitpath(im.dtfiles{1});
end


m = 1;
a = strrep(ns,'nsines.DP','grating.SFM');
if exist(a,'file') im.sffiles{n} = a;  m = m+1; end
a = strrep(ns,'nsines.DP','grating.SF');
if exist(a,'file') im.sffiles{n} = a;  m = m+1; end

if m > 1
    im.sffile = splitpath(im.sffiles{1});
end


function sta = CalcSTA(im, ids, xv)


if isfield(im,'dp')
    lp = im.phases - im.dp/2;
    rp = im.phases + im.dp/2;
    coss = im.sfs .* cos(lp);
    sins = sin(lp) .* im.sfs;
    if isfield(im,'rsfs')
        rcoss = im.rsfs .* cos(rp);
        rsins = sin(rp) .* im.rsfs;
    else
        rcoss = im.sfs .* cos(rp);
        rsins = sin(rp) .* im.sfs;
    end
else
    coss = im.sfs .* cos(im.phases);
    sins = sin(im.phases) .* im.sfs;
end
c = mean(coss(ids,:));
s = mean(sins(ids,:));
lsta = sc2space(im.sfvals, xv,c,s);
c = mean(rcoss(ids,:));
s = mean(rsins(ids,:));
rsta = sc2space(im.sfvals, xv,c,s);
sta = [lsta rsta];

function SetColorBar(cb)

if ~isempty(cb)
    pos = get(cb,'Position');
    pos(1) = pos(1) + 0.02;
    set(cb,'Position',pos,'YAxislocation','right');
end


function x = sc2space(sfs, deg, c,s)


for j = 1:length(sfs)
   cmps(j,:) =  c(j) .* cos(deg * sfs(j) * 2 * pi) + s(j) .* sin(deg * sfs(j) * 2 * pi);
end
x = mean(cmps);