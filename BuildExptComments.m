function [Comments, ExptList] = BuildExptComments(path, varargin)
%BuildExptComments(path, varargin) Extracts comment lines from
% Spike2.mat files and puts them in path/ExptComments.mat
% By default does nothing in ExptComments.mat exists.
rebuild = 0;
showcomments = 0;
j = 1;
while j <=length(varargin)
    if strncmpi(varargin{j},'rebuild',4)
        rebuild = 1;
    elseif strncmpi(varargin{j},'show',4)
        showcomments = 1;
    end
    j = j+1;
end

if ~isdir(path) && exist(path,'file')
    [path, filename] = fileparts(path);
end

outname = [path '/ExptComments.mat'];
if exist(outname) && ~rebuild
    load(outname);
else

d = mydir([path '/*idx.mat']);
nc = 1;
ne = 1;
Comments = [];
ExptsList = [];
for j = 1:length(d)
    load(d(j).name);
    for k = 1:length(Expt.Comments.text)
        if strncmp(Expt.Comments.text{k},' Electrode',7)
        elseif strncmp(Expt.Comments.text{k},'Electrode',7)
        elseif strncmp(Expt.Comments.text{k},'rf=',3)
        elseif strncmp(Expt.Comments.text{k},'cm=back=',8)
        else
            Comments(nc).text = Expt.Comments.text{k};
            Comments(nc).time = ConvertTime(Expt, Expt.Comments.times(k));
            Comments(nc).depth = FindElectrodeDepth(Expt, Expt.Comments.times(k));
            nc = nc+1;
        end
    end
    for k = 1:length(ExptList)
        ExptsList(ne).name = ExptList(k).expname;
        ExptsList(ne).start =  ConvertTime(Expt, ExptList(k).start./10000);
        ExptsList(ne).end =  ConvertTime(Expt, ExptList(k).end./10000);
        ExptsList(ne).depth(1) =  FindElectrodeDepth(Expt, ExptList(k).start);
        ExptsList(ne).depth(2) =  FindElectrodeDepth(Expt, ExptList(k).end);
        ExptsList(ne).ntrials =  CountTrials(Expt, [ExptList(k).start ExptList(k).end]);
        ne = ne+1;
    end
end
if ~isempty(Comments)
[a, b] = sort([Comments.time]);
Comments = Comments(b);
end
ExptList = ExptsList;
if ~isempty(ExptList)
[a, b] = sort([ExptList.start]);
ExptList = ExptList(b);
end
save(outname,'Comments','ExptList');
end


function truet = ConvertTime(Expt, t)
%t in seconds

truet = Expt.Header.CreationDate + t ./(24 .* 60 .* 60);

function d = FindElectrodeDepth(Expt, t)

bid = find(Expt.Trials.Start < t);
aid = find(Expt.Trials.Start > t);
if isempty(aid) && isempty(bid)
    id = 1;
elseif isempty(aid)
    id = bid(end);
elseif isempty(bid)
    id = aid(1);
else
    id = bid(end);
end
d = Expt.Trials.ed(id);

function nt = CountTrials(Expt, t)
bid = find(Expt.Trials.Start > t(1) & Expt.Trials.Start < t(2));
nt = sum(Expt.Trials.Result(bid) > 0);


