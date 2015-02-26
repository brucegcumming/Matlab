function DATA = cListExpts(DATA, Expts, varargin)

SpkDefs;
showwarn = 1;
explist = {};
na = 1;
nb = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nowarn',6)
        showwarn = 0;
    end
    j = j+1;
end
explist{1} = 'All';
explabel{1} = 'All';
exptnos = [];
if DATA.appending
    exptypelist = DATA.exptypelist;
    explist = DATA.explist;
    explabel = DATA.explabel;
    na = length(exptypelist)+1
else
    exptypelist = [];
    na= 2;
end
for j = 1:length(Expts)
    % expname is set in APlaySpkFile. Expt2Name n use unless it is made a
    % separate function to be used by both.
    %BUT exptype name IOS calculated fo real here, so they need to match!!
    %explist needs to mathc with what in Header.expname
    [expname, exptypename, suff] = Expt2Name(Expts{j});
    expname = [expname suff];
    exptypename = [exptypename suff];
    % if this expt is not in the list, add it
    eid = strmatch(Expts{j}.Header.expname,explist,'exact');
    if isempty(eid)
        explist{na} = expname;
        if DATA.state.online == 2
            explist{na} = Expts{j}.Header.expname;
            %            Expts{j}.Header.expname = expname;
        else
            explist{na} = Expts{j}.Header.expname;
        end
        exptypelist{na} = exptypename;
        DATA.explist = explist;
        DATA.exptypelist = exptypelist;
        cls = '';
        for k = 1:3
            outname = cmb.CombinedName(DATA, na,k);
            if exist(outname,'file')
                cls = [cls ' c' num2str(k)];
            end
        end
        if na > length(DATA.allcombineids)
            DATA.allcombineids{na} = [];
        end
        explabel{na} = [explist{na} cls];
        lens(na) = length(explist{na});
        if strcmp(expname,'cylinder.TwoCylDispXhxXbhP')
            explist{na} = [explist{na}];
            explist{na+1} = [explist{na}];
            explabel{na} = [explist{na} 'TWO' cls];
            explabel{na+1} = [explist{na+1} 'DID' cls];
            exptypelist{na} = [exptypename 'TWO'];
            exptypelist{na+1} = [exptypename 'DID'];
            na = na+1;
            DATA.exptypelist = exptypelist;
        end
        DATA.explist = explist;
        exids(j) = na;
        na = na+1;
    else
        exids(j) = eid;
        %        explist{na} = 'unknown';
    end
    DATA.exptnos(j) = GetExptNumber(Expts{j});
end


if isfield(DATA,'CellList')
    if ~isfield(DATA,'CellDetails')
        id = DATA.exptnos(find(exids ==j));
    else
        if size(DATA.CellList,1) < j
            str = [];
            if isfield(DATA,'CellDetails')
                eid = setdiff(DATA.exptnos,DATA.CellDetails.exptids);
                if ~isempty(eid);
                    str = sprintf('\nEx %s not in CellList',sprintf('%d ',eid));
                end
                if isfield(DATA,'suffixlist')
                    eid = setdiff(DATA.CellDetails.exptids,DATA.suffixlist);
                end
                if ~isempty(eid);
                    str = [str sprintf('\nEx %s not in Combine',sprintf('%d ',eid))];
                end
            end
            if showwarn
                acknowledge(['Number of Expts in CellList does not match' str],DATA.toplevel,'label','List Mismatch');
            end
        end
        for j = 2:length(explist)
            id = DATA.exptnos(find(exids ==j));
            id = find(ismember(DATA.CellDetails.exptids,id));
            ec = DATA.CellList(id,:,:);
            ec = ec(:);
            cellid = unique(ec(ec > 0));
            nspc = 2+round((max(lens)-length(explist{j}))*2);
            cls = [repmat(' ',1,nspc) 'Cells' sprintf(' %d',cellid)];
            explabel{j} = [explist{j} cls];
        end
    end
end

if isfield(DATA,'clst')
    p  = get(DATA.clst,'Listboxtop');
    if p > length(explabel)
        set(DATA.clst,'ListboxTop',1);
    end
    p  = get(DATA.clst,'value');
    if p > length(explabel)
        set(DATA.clst,'value',1);
    end
    set(DATA.clst,'string',explabel);
end
for j = 1:length(explist)
    DATA.combinenames{j} = cmb.CombinedName(DATA,j,0);
end
DATA.explist = explist;
DATA.exptypelist = exptypelist;
DATA.explabel = explabel;

