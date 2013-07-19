function FixSpkDir(name, varargin)

deletefile = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'delete',3)
        deletefile = 1;
    end
    j = j+1;
end

d = dir([name '/*.mat']);
mnk = GetMonkeyName(name);

mk = strmatch(mnk,{d.name});
if strfind(name,'/G')
oldf = strmatch('G',{d.name});
else
oldf = strmatch('M',{d.name});
end

for j = 1:length(mk)
    of = strrep(d(mk(j)).name,mnk,'');
    id = strmatch(of,{d(oldf).name});
    if length(id) == 1
        fprintf('%s matches %s\n',d(mk(j)).name,d(oldf(id)).name);
        if d(oldf(id)).datenum < d(mk(j)).datenum
            dname =[name '/' d(oldf(id)).name];
            fprintf('deleting %s\n',dname);
            if deletefile
                delete(dname);
            end
        elseif d(oldf(id)).datenum > d(mk(j)).datenum
            dname =[name '/' d(mk(j)).name];
            fprintf('deleting %s\n',dname);
            %delete(dname);
        end
    end
end