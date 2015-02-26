function str = Code2String(code, type)
%Code2String(code, type) return error string from codes
%Code2String(code, 'unsafetosave')  describes meaning of this field for a
%              reclassisy in AllVPcs

if nargin ==1
    type = 'unsafe';
end
if strncmpi(type,'unsafetosave',5)
    codes = {'Different Numbers of Classified Spikes ' ...
        'Correlation between classification space values low ' ...
        'Missing >10% of triggers ' };
    if code < -100 
        str = 'Cannot check - saved cluster not quantified';
        return;
    end
    
    if code < 0
        str = 'AutoCut ';
    else
        str = '';
    end
    code = abs(code);
    if bitand(code,1)
        str = [str codes{1}];
    end
    if bitand(code,2)
        str = [str codes{2}];
    end
    if bitand(code,4)
        str = [str codes{3}];
    end
end