function SaveConfig(DATA, file, savefields, varargin)
%SaveConfig(DATA, file, savefields, varargin)    
%saves state of fields in DATA (usually a gui UserData) to a file
%savefields is a cell array of fields that are saved
%Currently drops down one level only. ie.
% will save DATA.flag.next but not DATA.flag.next.next

verbose = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose',5)
        verbose = 1;
    end
    j= j+1;
end

fid = fopen(file,'w');
    if fid < 0
        fprintf('Cant write to %s\n',file)
        return;
    end

    for f = savefields;
        ndot = strfind(f{1},'\.')
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
                elseif isnumeric(val)
                    fprintf(fid,'DATA.%s.%s=[%s];\n',f{1},sf{j},num2str(val));
                elseif ischar(val)
                    fprintf(fid,'DATA.%s.%s=''%s'';\n',f{1},sf{j},val);
                end
            end
        end
    end
    if verbose
        fprintf('Saved Config to %s\n',file);
    end

