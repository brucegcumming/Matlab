function CloseAllHidden()

hh = allchild(0);
h = get(0,'Children');
    close(setdiff(hh,h));