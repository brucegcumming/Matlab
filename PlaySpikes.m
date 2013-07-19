function PlaySpikes(Spikes, varargin)

plotspikes = 0;
scale = 0;
j = 1;
if ischar(Spikes)
   set(findobj('Tag','SpikePause'),'value',2);
    return;
end
while j <= nargin - 1
    if strncmpi(varargin{j},'gui',3)
        gui= MakeGui;
        drawnow;
    elseif strncmpi(varargin{j},'plot',4)
        plotspikes = 1;
    end
        
    j = j+1;
end

if scale == 0
    scale = 1/max(max(abs(Spikes(2:33,:))));
end
lastspk = 1;
t = 1;
duration = 20000 * 3.2;
spk = 1;
paused = 0;
while spk < size(Spikes,2)
    tic;
    sound = zeros(1,duration);
    start = Spikes(1,lastspk);
    t = 1;
    if gui
        onesp = get(findobj('Tag','NextSpike'),'value');
        onetrial = get(findobj('Tag','NextTrial'),'value');
        paused = get(findobj('Tag','SpikePause'),'value');
        loop = get(findobj('Tag','SpikeLoop'),'value');
        if onesp
            paused = 2;
            set(findobj('Tag','NextSpike'),'value',0);
            drawnow;
        end
        if onetrial
            paused = 3;
            set(findobj('Tag','NextTrial'),'value',0);            
        end
        if isempty(onesp) %% window closed
            return;
        end
%%        fprintf('Paused %d %d\n',onesp,paused);
    end
    if paused == 1
        onesp = 0;
        pause(0.1);
    elseif paused == 2
        sound = zeros(1,3200);
        t = 1000;
        sound(t:t+31) = Spikes(2:33,spk);
        sound = sound * scale;
        wavplay(sound,32000);
        plot(Spikes(2:33,spk));
        hld = get(findobj('Tag','HoldSpike'),'value');
        if hld
            hold on;
        else
            hold off;
        end
        
        spk = spk+1;
        lastspk = spk;
        paused = 1;
    else
        while t < duration - 32 & spk < size(Spikes,2)
            sound(t:t+31) = Spikes(2:33,spk);
            spk = spk + 1;
            t = ceil( (Spikes(1,spk) - start) * 3.2);
        end
        if loop & spk >= size(Spikes,2)
            spk = 1;
            lastspk = 1;
        else
            lastspk = spk;
        end
    toc;
    sound = sound * scale;
    if strcmp(computer,'PCWIN')
        if plotspikes
            plot(sound);
            drawnow;
        end
        wavplay(sound,32000);
        if gui
            set(gui, 'Name',sprintf('Spike %d',spk));
        end
        pause(0.1);
        if paused == 3
            paused = 1;
        end
    else
        plot(sound);
    end
    end
end


function restart
spk = 1;

function cntrl_box = MakeGui(varargin)

cw=10;
SPACE = 10;
VSPACE = 10;
cntrl_box = figure('Position', [100 100  200 200],...
    'NumberTitle', 'off', 'Tag','SpikeGui','Name','Spike Playback');
  bp(1) = SPACE; bp(2) = VSPACE; bp(3) = 80; bp(4) = 22;

    hm = uimenu(cntrl_box,'Label','Mark');
%  uimenu(hm,'Label','&Restart','Callback',&restart);

bp(2) = bp(2)+bp(4);
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Pause', 'Tag', 'SpikePause', 'Position', bp,'value',1);

bp(2) = bp(2)+bp(4);
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Next', 'Tag', 'NextSpike', 'Position', bp);

bp(2) = bp(2)+bp(4);
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Next Trial', 'Tag', 'NextTrial', 'Position', bp);

bp(2) = bp(2)+bp(4);
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Loop', 'Tag', 'SpikeLoop', 'Position', bp);

bp(2) = bp(2)+bp(4);
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Overlay', 'Tag', 'HoldSpike', 'Position', bp);
