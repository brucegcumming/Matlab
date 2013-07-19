function Expt = CheckEyeData(Expt)
%Expt = CheckEyeData(Expt) makes sure lengths of EyeData
%in Trials match. 


for j = 1:length(Expt.Trials)
    emlen(j) = length(Expt.Trials(j).EyeData);
    if emlen(j) > 0
        if isnan(Expt.Trials(j).emstarti)
            prepts(j) = 0;
        else
            prepts(j) = Expt.Trials(j).emstarti;
        end
        postpts(j) = emlen(j)-prepts(j);
    else
        prepts(j) = 0;
        postpts(j) = 0;
    end
end

pre = max(prepts)
post = max(postpts);
for j = 1:length(Expt.Trials)
    prepad = zeros(pre-prepts(j),4);
    postpad = zeros(post-postpts(j),4);
    E = cat(1, prepad, Expt.Trials(j).EyeData, postpad);
    Expt.Trials(j).EyeData = E;
end
emlen = post+pre;
Expt.Header.emlen = post+pre;
Expt.Header.emtimes = (([1:emlen]-pre) .* 1./Expt.Header.CRsamplerate);

min(emlen);