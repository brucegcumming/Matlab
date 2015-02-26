function s = nonewline(s)

while  ~isempty(s) && s(end) == 10
    s = s(1:end-1);
end