function TextText(name)

method = 1;

load(name);
Text = Ch30;

if sum(method == 1)
tic;
id = strmatch(Text.text,'xxx=');
ids = setdiff(1:size(Text.text,1),id);
xid = ones(size(id));
nl = size(Text.text,1)-length(id);
aText.text = mat2cell(Text.text(ids,1:end-1),ones(1,length(ids)),size(Text.text,2)-1)
k=1;
m = 1;
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
        xid(m) = j;
        m = m+1;
        if strncmp(aText.text{k},'mtcL',4)
            m 
        end
% if we ever go back to this, probably need deblank(Text.....
    end
    if length(aText.text{k}) == 0 
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
    end
    k = k+1;
end
aText.times = Text.times(ids);
aText.codes = Text.codes(ids,:);
aText.text = deblank(aText.text);
toc;
end

if sum(method == 2)
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        xid(k) = j;
        k = k+1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
    else
% if we ever go back to this, probably need deblank(Text.....
    aText.text{k} = Text.text(j,1:end-1);
    aText.times(k) = Text.times(j);
    aText.codes(k,:) = Text.codes(j,:);
    end
    if length(aText.text{k}) %%? safe may remove blank strings that have codes.
        k = k+1;
    else
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
        k = k+1;
    end
end
aText.text = deblank(aText.text);
fprintf('Reading Text %.2f\n',toc);
end

