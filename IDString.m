function s = IDString(X, varargin)

s = [];
type = [];
j = 1; 
while j <= length(varargin)
    if strncmp(varargin{j},'probe',4)
        type = probe;
    end
    j = j+1;
end

if isfield(X,'Header') && (isfield(X.Header,'cellnumber') || isfield(X.Header,'probe')) 
    P = X(1).Header; %in case its an array
    if ~isfield(P,'cellnumber')
        P.cellnumber = 0;
    end
    if isempty(type) || strcmp(type,'probe')
        if isfield(P,'probe')
            pstr = sprintf(' P%.0f',P.probe);
            if mod(P.probe,1) > 0.01
                pstr = sprintf('P%.1f',P.probe);
            end
        elseif isfield(X,'probe')
            pstr = sprintf(' P%.0f',X.probe);
        else
            pstr = '';
        end
        if P.cellnumber > 0
            s = sprintf('Cell %d%s',P.cellnumber, pstr);
        else
            s = sprintf('MU%s', pstr);
        end
    end
end