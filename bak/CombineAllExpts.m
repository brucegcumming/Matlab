function CombineAllExpts(dira, dirb, varargin)

da = dir(dira);
db = dir(dirb);
verbose = 0;

for j = 1:length(da)
    id = regexp(da(j).name,'\..*\.[a-z]*\.[a-z]*[A-Z]*\.mat');
    if ~isempty(id)
        suffix = da(j).name(id(1):end);
        if verbose
            fprintf('Expt: %s Suffix %s\n',da(j).name,suffix);
        end
        for k = 1:length(db)
            if strfind(db(k).name,suffix)
                outname = strrep(da(j).name,suffix,['C' suffix]);
                fprintf('Combine: %s %s -> %s\n',da(j).name,db(k).name,outname);
                a = load([dira '/' da(j).name]);
                b = load([dirb '/' db(k).name]);
                if isfield(a,'LFP') & isfield(b,'LFP')
                    LFP = CombineExpts(a.LFP,b.LFP);
                    save([dira '/' outname],'LFP');
                elseif ~isfield(a,'Expt')
                    fprintf('%s has no Expt\n',da(j).name);
                elseif ~isfield(b,'Expt')
                    fprintf('%s has no Expt\n',db(k).name);
                else
                Expt = CombineExpts(a.Expt,b.Expt);
                save([dira '/' outname],'Expt');
                end
            end
        end

end

end
        