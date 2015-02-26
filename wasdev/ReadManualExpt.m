function ReadManualExpt(name,varargin)


prefix = [name '/stim'];
stim=0;
txt = scanlines(sprintf('%s%d',prefix,stim));
while ~isempty(txt)
    fprintf('%d:',stim);
    for j = 1:length(txt)
        fprintf('%s ',txt{j});
    end
    fprintf('\n');
    stim = stim+1;
    txt = scanlines(sprintf('%s%d',prefix,stim),'silent');
end
    