function plotdata = StoreData(OTTF, data, id)

plotdata = [];
if isempty(data) & isempty(id)
    plotdata = get(findobj('Tag',OTTF.tag.data),'UserData');    
elseif id == 0
    set(findobj('Tag',OTTF.tag.data),'UserData',data);
elseif OTTF.state.storedata
    plotdata = get(findobj('Tag',OTTF.tag.data),'UserData');
    plotdata{id} = data;
    set(findobj('Tag',OTTF.tag.data),'UserData',plotdata);
end