function FindFrameTimes(Trials, frametimes)

ts = now;
if ~isfield(Trials,'End')
    Trials.End = Trials.Start + 20000;
end

ff = 1;
nf = length(frametimes);
maxl = 100+ceil(diff(Trials.Start)./100); 
maxl(end+1) = 100;
for j = 1:length(Trials.Start)
      if ~isempty(frametimes) % have VTR channel
%                id = find(frametimes > Trials.Start(j) & frametimes < Trials.Start(j)+500);
          if ff+maxl(j) > nf
              last = nf;
          else
              last = ff+maxl(j);
          end
          id = find(frametimes(ff:last) > Trials.Start(j),1);
          if ~isempty(id)
              ff = id(1);
          end
          
          id = find(frametimes(ff:ff+1000)> Trials.End(j),1);
      end
end
    
fprintf('Took %.2f\n',mytoc(ts));