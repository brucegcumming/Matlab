function I = BuildLFPIndex(name,varargin)
%I = BuildLFPIndex(name,varargin)
%For Utah Folders, read in all the lfp files that are made and summarize contents
%For Spike2 Folders, compare Expts List with list of .lfp files
listonly = 0;
j = 1;
defaultarray = '';

while j <= length(varargin)
    if strncmpi(varargin{j},'fixlist',4)
        listonly = 2;
    elseif strncmpi(varargin{j},'list',4)
        listonly = 1;
    end
    j = j+1;
end

if iscellstr(name)
    for j = 1:length(name)
        I{j} = BuildLFPIndex(name{j},varargin(:));
    end
    return;
end


if isdir(name)
    mnk = GetMonkeyName(name);
    if isempty(defaultarray)
        if strcmp(mnk,'jbe')
            defaultarray = 'Utah';
        else
            defaultarray = 'uprobe';
        end
    end
    goodexpts = 0;
    outname = [name '/LFPlist.mat'];
    Array = GetArrayConfig(name);
    if isempty(Array)
            Array.type = defaultarray;
    end
    if listonly 
        if sum(strncmpi(Array.type,{'uprobe' 'vprobe'},4))
            Expts = ReadExptDir(name);
            for j = 1:length(Expts)
                lfpname = [name '/' GetName(Expts{j},'lfp')];
                if ~exist(name)
                    fprintf('Missing %s\n',name);
                    goodexpts(j) = 0;
                else
                    goodexpts(j) = 1;
                end
            end
            I.name = outname;
            I.havefiles = goodexpts;
            I.missing = sum(goodexpts ==0);
            I.ngood = sum(goodexpts > 0);
        elseif strncmpi(Array.type,'utah',4)
            I.name = outname;
            if ~exist(outname)
                I.list = 0;
                fprintf('No LFP List for %s\n',name);
            else
                I.list = 1;
                I.exptno = [];
                I.lowpasserror = 0;
                OUT = load(outname);
                if ~iscell(OUT.D)
                    D{1} = OUT.D;
                    OUT.D = D;
                end
                for j = 1:length(OUT.D)
                    exptno = 0;
                    ncut = 0;
                    
                    if isfield(OUT.D{j},'exptno')
                        name = [outname 'Expt' num2str(OUT.D{j}.exptno)];
                        exptno = OUT.D{j}.exptno;
                    end
                    I.exptno = [I.exptno exptno];
                    if isfield(OUT.D{j},'lfpnames')
                        I.names{j} = OUT.D{j}.lfpnames{1};
                        name = OUT.D{j}.lfpnames{1};
                    elseif exptno > 0
                        name = [outname 'Expt' num2str(exptno)];
                    else
                        name = outname;
                    end
                    if isfield(OUT.D{j},'errs')
                        I.nerrs(j) = length(OUT.D{j}.errs);
                        for e = 1:length(OUT.D{j}.errs)
                            err = OUT.D{j}.errs{e};
                            if strfind(err,'points, array is')
                                fprintf('Bad Channel in %s\n',name);
                            elseif strfind(err,'filtered out')
                                ncut = ncut+1;
                            else
                                fprintf('Unknown Errors in Expt %d\n',j);
                            end
                        end
                    else
                        I.nerrs(j) = 0;
                    end
                    I.lowpasserror(j) = ncut;
                end
                if  sum(I.lowpasserror)
                    fprintf('Filteriing wrong in %s\n',name);
                end
            end
        else
            I.err = sprintf('Unknown Array Type %s',Arrat,type);
            cprintf('red','%s in %s\n',I.err,name);
            I.name = outname;
        end
        if listonly == 1
            return;
        end
    else
        I.name = outname;
    end
    if exist(outname)
        OUT = load(outname);
        if isfield(OUT,'D')
            if  isstruct(OUT.D)
                olddata = OUT.D;
            OUT.D = {};
            else
                D = {};
                for j = 1:length(OUT.D)
                    expts(j) = ExptIndex(OUT.D{j});
                    if expts(j) > 0
                        D(expts(j)) = OUT.D(j);
                    end
                end
                OUT.D = D;
            end
        end
        
    end
    if sum(strncmpi(Array.type,{'uprobe' 'vprobe'},4))
        Expts = ReadExptDir(name);
        I.LFPerrs = {};
        I.LFPerrdata =[];
        for j = 1:length(Expts)
            [~, E] = LoadLFP(Expts{j});
            if E.gotlfp == 0
                goodexpts(j) = 0;
            elseif isfield(E,'LFPerrs')
                goodexpts(j) = 2; %have data, but with errors
                I.LFPerrs = cat(1,I.LFPerrs(:),E.LFPerrs(:));
                I.LFPerrdata = cat(2,I.LFPerrdata,E.LFPerrdata);
            else
                goodexpts(j) = 1;
            end
        end
        I.goodexpts = goodexpts;
        if sum(goodexpts)
            save(outname,'-struct','I');
        end
    elseif strncmpi(Array.type,'Utah',4) 
        d = mydir([name '/Expt*.lfp.mat']);
        for j = 1:length(d)
            X = load(d(j).name);
            e = GetExptNumber(d(j).name);
            if isfield(X,'LFPHeader')
                OUT.D{e}.Header = X.LFPHeader;
                goodexpts(e) = 1;
            elseif isfield(X,'LFP')
                OUT.D{e}.npts = size(X.LFP.rawlfp,1);
                if isfield(X.LFP,'Header')
                    OUT.D{e}.Header = X.LFP.Header;
                else
                    OUT.D{e} = CopyFields(OUT.D{e},X.LFP,'kernel','sd','decimate','samper');
                end
                goodexpts(e) = 1;
            end
            OUT.D{e}.exptno = e;
        end
        if sum(goodexpts)
            bid = find(~goodexpts); %no data here any more
            for j =1:length(bid)
                OUT.D{bid(j)} = [];
            end
            save(outname,'-struct','OUT');
        end
    else
        cprintf('red','Unknown Array Type %s in %s\n',Array.type,name);
    end
else
    I.name = name;
end

function eid = ExptIndex(X)

eid = 0;
if isfield(X,'exptno')
    eid = X.exptno;
elseif isfield(X,'Header')
end
