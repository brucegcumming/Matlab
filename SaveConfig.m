function file = SaveConfig(DATA, file, savefields, varargin)
%SaveConfig(DATA, file, savefields, varargin)    
%saves state of fields in DATA (usually a gui UserData) to a file
%savefields is a cell array of fields that are saved
%Currently drops down one level only. ie.
% will save DATA.flag.next but not DATA.flag.next.next

verbose = 0;
interactive = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose',5)
        verbose = 1;
    elseif strncmpi(varargin{j},'choose',5)
        interactive = 1;
    end
    j= j+1;
end

if interactive
    [name, outdir] = uiputfile(file);
    if ~ischar(name)
        return;
    end
    file = [outdir '/' name];
end
BackupFile(file);
fid = fopen(file,'w');
    if fid < 0
        cprintf('errors','Cant write to %s\n',file)
        return;
    end

    for f = savefields;
        ndot = strfind(f{1},'\.');
        if isempty(strfind(f{1},'\.'))
            if isfield(DATA,f{1})
            val = DATA.(f{1});
            end
        elseif length(ndot) == 1
            a = f{1}(1:ndot-1);
            b = f{1}(ndot+1:end);
            val = DATA.(a).(b);
        end
        if isempty(val)
            fprintf(fid,'DATA.%s=[];\n',f{1});
        elseif isnumeric(val)
            fprintf(fid,'DATA.%s=[%s];\n',f{1},num2str(val));
        elseif ischar(val)
            fprintf(fid,'DATA.%s=''%s'';\n',f{1},val);
        elseif isstruct(val)
            sf = fields(val);
            for j = 1:length(sf)
                val = DATA.(f{1}).(sf{j});
                if isempty(val)
                    fprintf(fid,'DATA.%s.%s=[];\n',f{1},sf{j});
                elseif isnumeric(val) || islogical(val)
                    fprintf(fid,'DATA.%s.%s=[%s];\n',f{1},sf{j},num2str(val));
                elseif ischar(val)
                    fprintf(fid,'DATA.%s.%s=''%s'';\n',f{1},sf{j},val);
                end
            end
        elseif iscell(val)
            for j = 1:size(val,1);
                for k = 1:size(val,2)
                    cval = DATA.(f{1}){j,k};
                    if ischar(cval)
                        fprintf(fid,'DATA.%s{%d,%d}=''%s'';\n',f{1},j,k,cval);
                    elseif isstruct(cval)
                        for p = fields(cval)'
                            if ishandle(cval.(p{1}))
                                fprintf(fid,'DATA.%s{%d,%d}.%s=%f;\n',f{1},j,k,p{1},double(cval.(p{1})));
                            else
                                fprintf(fid,'DATA.%s{%d,%d}.%s=%f;\n',f{1},j,k,p{1},cval.(p{1}));
                            end
                        end
                    else
                        fprintf(fid,'DATA.%s{%d,%d}=%f\n',f{1},j,k,cval);
                    end
                end
            end
        end
    end
    if verbose
        fprintf('Saved Config to %s\n',file);
    end

fclose(fid);