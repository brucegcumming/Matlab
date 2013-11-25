function X = myNormalize(X, varargin)
%Y = myNormalize(X,...)
%scale each row of a Matrix
%default is to divide each row by max(row);
mode = 'max';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'mean',3)
        mode = 'mean';
    elseif strncmpi(varargin{j},'std',3)
        mode = 'std';
    end
    j = j+1;
end

for j = 1:size(X,1)
    if strcmp(mode,'std')
        m = std(X(j,:));
    elseif strcmp(mode,'mean')
        m = mean(X(j,:));
    else
        m = max(X(j,:));
    end
        X(j,:) = X(j,:)./m;
end
