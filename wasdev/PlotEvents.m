function PlotEvents(t, pre, post, varargin)
Events.times = [];
StimOn.times = [];
frames.times = [];
SpkDefs;

for j = 1:length(varargin)
   if isstruct(varargin{j}) & isfield(varargin{j},'times')
       if isfield(varargin{j},'text')
           Text = varargin{j};
       elseif strncmp(varargin{j}.title,'VTR',3)
           frames = varargin{j};
       elseif strncmp(varargin{j}.title,'StimOn',6)
           StimOn = varargin{j};
       elseif isfield(varargin{j},'codes')
           Events = varargin{j};
       end
   end
end


hold off; 
for j = 1:length(t)
    line = j-1;
    ts = [t(j) - pre t(j) + post];
    id = find(StimOn.times > ts(1) & StimOn.times < ts(2));
    if length(id)
    ltimes = sort(cat(1,StimOn.times(id),StimOn.times(id)));
    for k = 1:2:length(ltimes)
        level(k) = StimOn.level(ceil(k/2));
        level(k+1) = ~(StimOn.level(ceil(k/2)));
    end
    plot(ltimes,level,'r');
    hold on;
    end
    id = find(Events.times > ts(1) & Events.times < ts(2) & Events.codes(:,1) == FRAMESIGNAL);
    plot(Events.times(id),ones(size(id)).*0.9,'b+','markersize',3);
    hold on;
    id = find(Events.times > ts(1) & Events.times < ts(2) & Events.codes(:,1) == ENDSTIM);
    plot(Events.times(id),ones(size(id)).*0.8,'m+','markersize',3);
    plot([t(j) t(j)],[0.4 0.6]+line,'k','linewidth',2)
end

       