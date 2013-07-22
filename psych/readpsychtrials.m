function [Trials, ExptTrialList] = readpsychtrials(file, varargin)

j = 1;
k = 1;
doall = 0;
expctr = 1;
savefile = 0;
xo = [];
ypos = [];
varargon = {};
while(j < nargin)
    if(strncmpi(varargin{j},'all',3))
        doall = 1;
    elseif(strncmpi(varargin{j},'save',4))
        savefile = 1;
    else
        varargon{k} = varargin{j};
        k = k+1;
    end
    j = j+1;
end

expstarts = [];
expno = 1;
probit = [];
or = [];
ndat = 1;
fid = fopen(file);
tline = fgetl(fid);
nt = 1;
sortcode = 'sd';
sd = 0;
sign = 0;
while(ischar(tline))
    if(length(tline) > 0)
  if(tline(1) == 'R' & ismember(tline(2),'WGB'))
      [c, count, err, cpos]  = sscanf(tline,'R%c ',1);
      if(c == 'G')
          Trials(nt).score = 1;
      elseif(c == 'W')
          Trials(nt).score = -1;
      elseif(c == 'B')
          Trials(nt).score = 0;
      end
      id = strfind(tline,'sn=');
      if ~isempty(id)
          sign = sscanf(tline(id:end),'sn=%f');
          
      end
      [Trials(nt).x, count, err, xpos]  = sscanf(tline(cpos+3:end),'%f ',1);
      if(~strncmp(tline(cpos+xpos+2:end),'fl=',3))
          [Trials(nt).x2, count, err, xpos]  = sscanf(tline(cpos+xpos+5:end),'%f ',1);
      end
      if(length(sortcode))
          Trials(nt).sd = sd;
      end
      Trials(nt).xo = xo;
      Trials(nt).or = or;
      Trials(nt).yo = ypos;
      Trials(nt).id = sid;
      if(Trials(nt).score ~= 0 | doall)
          Trials(nt).sign = sign;
          nt = nt+1;
      end
  elseif(strncmp(tline,'id',2))
      sid = sscanf(tline,'id%f');
  elseif(strncmp(tline,'xo',2))
      xo = sscanf(tline,'xo%f');
  elseif(strncmp(tline,'yo',2))
      ypos = sscanf(tline,'yo%f');
  elseif(strncmp(tline,'sq',2))
      sign = sscanf(tline,'sq%f');
  elseif(strncmp(tline,'vs',2))
      sign = sscanf(tline,'vs%f');
  elseif(strncmp(tline,'sn=',3))
      sign = sscanf(tline,'sn=%f');
  elseif(strncmp(tline,'et',2) && ~strncmp(tline,'et ',3))
   if expctr < 2 | nt > expstarts(expctr-1)
       expstarts(expctr) = nt;
       exptype{expctr} = sscanf(tline,'et%2s');
       expctr = expctr+1;
   end
  elseif(strncmp(tline,'e2',2))
       expbtype{expctr} = sscanf(tline,'e2%2s');
  elseif(strncmp(tline,'Stimulus',4))
   if expctr < 2 | nt > expstarts(expctr-1)
    expstarts(expctr) = nt;
    expctr = expctr+1;
   end
   id = strfind(tline,'xo=');
   if ~isempty(id)
       xo = sscanf(tline(id:end),'xo=%f');
   end
   id = strfind(tline,'or=');
   if ~isempty(id)
       or = sscanf(tline(id:end),'or=%f');
   end
  elseif(strncmp(tline,sortcode,2)) 
    sd  = sscanf(tline,'%*c%*c%f');
  end
end
  tline = fgetl(fid);
end


ExptTrialList.n = expstarts;
ExptTrialList.exptypes = exptype;
ExptTrialList.expbtypes = expbtype;
ExptTrialList.blocks = expstarts;
for j = 1:length(expstarts)
    ExptTrialList.times{j} = sprintf('%d',expstarts(j));
    ExptTrialList.names{j} = sprintf('%d',expstarts(j));
end

if(savefile)
    save(file,'Trials','ExptTrialList');
end
