function X = myNormalize(X, varargin)
%Y = myNormalize(X,...)
%scale each row of a Matrix
%default is to divide each row by max(row);
%Y = myNormalize(X,S) where S is a vector
%           Divides X(j,:) by S(j);

mode = 'max';
scales = [];
j = 1;
while j <= length(varargin)
    if isnumeric(varargin{j}) && length(varargin{j}) == size(X,1) 
    elseif strncmpi(varargin{j},'mean',3)
        mode = 'mean';
    elseif strncmpi(varargin{j},'std',3)
        mode = 'std';
    end
    j = j+1;
end

for j = 1:size(X,1)
    if ~isempty(scales)
        X(j,:) = X(j,:)./scales(j);
    else
    if strcmp(mode,'std')
        m = std(X(j,:));
    elseif strcmp(mode,'mean')
        m = mean(X(j,:));
    elseif strcmp(mode,'user')
        m = mean(X(j,:));
    else
        m = max(X(j,:));
    end
        X(j,:) = X(j,:)./m;
end
end
