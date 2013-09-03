function [err, counts] = CheckExptRates(Expt, varargin)

printwarn = 0;
warncolor = 'red';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'print',5)
        printwarn = 1;
    end
    j = j+1;
end

    if iscell(Expt)
        [err, counts] = CheckExptsRateSequence(Expt, varargin{:});
        return;
    end
    if isfield(Expt,'Spikes') %An AllExpt Structure
        AllExpt = Expt;
        for j = 1:length(AllExpt.Spikes)
            [err{j}, counts{j}] = CheckExptRates(All2Expt(AllExpt,j),varargin{:});
        end
        return;
    end
    err.warning = zeros(1,5);
    counts = [Expt.Trials.count];
    if length(counts) > 100
        smw = 10;
    elseif length(counts) > 20
        smw = 5;
    elseif length(counts) == 0
        smw = 2;
    else
        smw = 2;
    end
    idstr = Expt2Name(Expt,'addprobe');

    smc = smooth(counts, smw);
    err.ff = var(smc)./mean(smc);
    b = polyfit(1:length(counts),sqrt(counts),1);
    err.mean = mean(counts);
    err.slope = length(counts).*b(1)./mean(sqrt(counts));
    err.smw = smw;
    if isfield(Expt.Header,'exptno')
        err.exptno = Expt.Header.exptno;
    end
    if sum(counts) == 0
        err.warning(6) = 1;
        mycprintf('blue','E%dCell%d no spikes (%s)\n',err.exptno,Expt.Header.cellnumber,idstr);
        err.slope = 0;
    end
    if err.ff > 5
        err.warning(4) = 1;
        if printwarn
            mycprintf(warncolor,'%s:Fano Factor %d pt avg %.1f\n',idstr,smw,err.ff);
        end
    end
    if isfield(Expt.Header,'BlockStart')
        trls = [Expt.Trials.Trial];
        for b = 1:length(Expt.Header.BlockStart)
            if b < length(Expt.Header.BlockStart)
                last = Expt.Header.BlockStart(b+1);
            else
                last = max(trls)+1;
            end
            id = find( trls >= Expt.Header.BlockStart(b) & trls < last);
            blkmean(b) = mean([Expt.Trials(id).count]);
            blkn(b) = length(id);
        end
        err.blkff = var(blkmean)./mean(blkmean);
        err.blkmean = blkmean;
        err.blkn = blkn;
        err.blkskew = skewness(blkmean);
        if err.blkff > 10
            err.warning(5) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Fano Factor %.1f\n',idstr,min(err.blkn),max(err.blkn),err.blkff);
            end
        end
        if err.blkskew < -1.5
            err.warning(2) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Skewness %.1f\n',idstr,min(err.blkn),max(err.blkn), err.blkskew);
            end
        end
    end
    
function [err, counts] = CheckExptsRateSequence(Expts,varargin)
printwarn = 0;
warncolor = 'red';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'print',5)
        printwarn = 1;
    end
    j = j+1;
end
err.warning = zeros(1,5);

    if iscell(Expts)
        nx = 0;
        for j = 1:length(Expts)
            if ~isempty(Expts{j})
                nx = nx+1;
                [ err.errs(nx), counts{nx}] = CheckExptRates(Expts{j});
                m(nx) = mean(counts{nx});
                err.blkn(nx) = length(counts{nx});
                idstr = Expt2Name(Expts{j},'addprobe');
            end
        end
        err.blkmean = m;
        id = find(m > 0);
        mrate = mean(cat(2,counts{id}));
        medrate = median(m(m>0));
        sd = std(cat(2,counts{id}));
        err.blkcv = std(m(id))./mean(m(id));
        err.blkff = var(m(id))./mean(m(id));
        err.blkskew = skewness(m); %N.B include all for this
        for j = 1:length(m)
            if m(j) > 0;
             err.diffs(j) = ( err.errs(j).mean)./medrate;
            else
             err.diffs(j) = 1;
            end
        end
        if  err.blkff > 5 || err.blkcv > 1
             err.warning(5) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Fano Factor %.1f CV %.2f\n',idstr,min(err.blkn),max(err.blkn), err.blkff,err.blkcv);
            end
        end
       
        if err.blkskew < -1.5
            err.warning(2) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Skewness %.1f\n',idstr,min(err.blkn),max(err.blkn), err.blkskew);
            end
        end
        if sum(m ==0) % should not happen for defined cell
            err.warning(3) = 1;
        end
    end

    
    
    
