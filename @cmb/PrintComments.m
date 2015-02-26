function DATA = PrintComments(DATA)

if isfield(DATA.Expts{1}.Header,'rfstr');
    fprintf('%s\n',DATA.Expts{1}.Header.rfstr);
end
if isfield(DATA,'Comments')
    if isfield(DATA.Comments,'Peninfo')
        fprintf('%s\n',DATA.Comments.Peninfo.trode);
        a= InterpretLine(DATA.Comments.Peninfo.trode);
        DATA.Header.trode = DATA.Comments.Peninfo.trode;
        CopyFields(DATA.Header,a);
    end
    for j = 1:length(DATA.Comments.text)
        if isempty(strmatch('cm=back=',DATA.Comments.text{j})) & ...
                isempty(strfind(DATA.Comments.text{j},'Tube Protrudes'))
            id = find(DATA.AllData.TrialTimes < DATA.Comments.times(j));
            if isempty(id)
                id = 1;
            end
            fprintf('%.1f(Trial%.0f): %s\n',DATA.Comments.times(j)./10000,...
                DATA.AllData.Trialids(id(end)),...
                DATA.Comments.text{j});
        end
    end
end
cmfile = strrep(DATA.datafilename,'.mat','.txt');
fid = fopen(cmfile,'r');
if fid > 0
    fprintf('Offline Comments:\n');
    a = textscan(fid,'%s','delimiter','\n');
    for j = 1:length(a{1})
        fprintf('%s\n',a{1}{j});
    end
    fclose(fid);
end


