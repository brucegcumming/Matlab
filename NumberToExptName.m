    function name = NumberToExptName(x, Expt)
        name = 'e0';
        switch x
            case 107
            case 112.5
            name = 'sf';
            case 101
            name = 'me';
            case 96
            name = 'Pp';
            case 106
            case 100
            name = 'Op';
            case 111
                name = 'or';
            otherwise
            fprintf('Unknown');
        end
        
        fprintf('Ex %d ->%s\n',x,name);