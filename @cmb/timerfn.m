function timerfn(tim, varargin)
DATA = get(findobj('Tag',get(tim,'Tag')),'UserData');
d = dir(DATA.datafilename);
if d.datenum > DATA.lastread || d.bytes > DATA.lastsize;
%       fprintf('Reading from %d\n',DATA.linesread);
[Expts, DATA.linesread] = ReadOnlineTxt(DATA.datafilename, DATA);
%        fprintf('Read %d expts\n',length(Expts));
Expts = cmb.CountTxtSpikes(Expts,DATA.probe,DATA.spikelist);
%         fprintf('List\n');
if length(Expts) > length(DATA.Expts)
DATA.Expts = Expts;
DATA = cmb.cListExpts(DATA,Expts);
%            fprintf('listexp\n');
DATA = cmb.combine('listexps',DATA,'Tag',DATA.tag.top);
%             fprintf('Set\n');
set(DATA.elst,'value',length(get(DATA.elst,'string')))
fprintf('Replot (%d,%d) at %s\n',length(Expts),DATA.linesread,datestr(now));
else
DATA.Expts = Expts;
fprintf('Replot (%d,%d,%d) at %s\n',length(Expts),length(Expts{end}.Trials),DATA.linesread,datestr(now));
end
DATA.lastread = d.datenum;
DATA.lastsize = d.bytes;
cmb.combine('setexp',DATA,'Tag',DATA.tag.top);
%        fprintf('Expt Set\n');
end


