function spikelist = WhichClusters(top,varargin)

spikelist = [];

if get(findobj(top,'Tag','UseCluster0'),'value')
    spikelist = [spikelist 0];
end
if get(findobj(top,'Tag','UseCluster1'),'value')
    spikelist = [spikelist 1];
end
if get(findobj(top,'Tag','UseCluster2'),'value')
    spikelist = [spikelist 2];
end
if get(findobj(top,'Tag','UseCluster3'),'value')
    spikelist = [spikelist 3];
end
if get(findobj(top,'Tag','UseCluster4'),'value')
    spikelist = [spikelist 4];
end
if isempty(spikelist)
    spikelist = -1;
end




