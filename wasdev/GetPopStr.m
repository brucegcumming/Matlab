function s = GetPopStr(a)
%getpopstr(x)  returns the curretly selected string value from a pop-menu
str = get(a,'string');
s = deblank(str(get(a,'value'),:));