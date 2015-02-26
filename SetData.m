function SetData(DATA) 
%SetData(DATA) shortcut for set(DATA.toplevel,'UserData',DATA);
interactive =0;
if isfield(DATA,'toplevel')
    set(DATA.toplevel,'UserData',DATA);
elseif interactive
    cprintf('error','No Field toplevel in DATA');
end