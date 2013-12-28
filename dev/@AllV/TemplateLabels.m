function Labels = TemplateLabels(DATA, usestd)    ispk = DATA.probe(1);    if isfield(DATA,'chspk')        chspk = DATA.chspk;    else        chspk = DATA.probe(1)+ [-1:1];        chspk = chspk(chspk >0 & chspk <= DATA.allnprobes);    end    if usestd    Labels{1} = sprintf('1r');    Labels{2} = sprintf('2r');    Labels{3} = sprintf('1dt');    Labels{4} = sprintf('2dt');    Labels{5} = sprintf('?');    Labels{6} = sprintf('?');    Labels{7} = sprintf('?');    Labels{8} = sprintf('2dt');    Labels{9} = sprintf('?');    Labels{10} = sprintf('2dt');    Labels{11} = sprintf('?');    Labels{12} = sprintf('?');    return;    end    if length(DATA.chspk) == 1        ispk = DATA.probelist(DATA.probe(1));        Labels{1} = sprintf('%d:r',ispk);        Labels{2} = sprintf('MU');        Labels{8} = sprintf('%d:dt',ispk);        Labels{9} = sprintf('Std1');        Labels{4} = sprintf('std2');        Labels{6} = sprintf('%d:dp',ispk);        Labels{5} = sprintf('%dr-otherr',ispk);        Labels{7} = sprintf('dp');        Labels{3} = sprintf('mu:dvdt');        Labels{10} = sprintf('dp wieghted');        Labels{11} = sprintf('trigger');        Labels{12} = sprintf('dvdt Std1',ispk);        Labels{14} = sprintf('MUsum',ispk);        Labels{15} = sprintf('sum2',ispk);        Labels{13} = sprintf('absdiff',ispk);        Labels{16} = sprintf('sum2-sum1');        Labels{17} = sprintf('sumP%d',DATA.probe(1)-1);        Labels{18} = sprintf('sumP%d',DATA.probe(1)+1);    else        Labels{1} = sprintf('%d:r',ispk);        Labels{2} = sprintf('sum');        Labels{8} = sprintf('%d:dt',ispk);        Labels{9} = sprintf('%d:dy',ispk);        Labels{5} = sprintf('%d:csd',chspk(1));        Labels{6} = sprintf('%d:dp',ispk);        Labels{7} = sprintf('%d:dp',chspk(1));        Labels{3} = sprintf('%d:r',chspk(1));        if length(chspk) > 2            Labels{4} = sprintf('%d:r',chspk(end));        else            if max(chspk)  == DATA.nprobes                xspk = min(chspk)-1;            else                xspk = max(chspk)+1;            end            Labels{4} = sprintf('%d:r',xspk);        end        Labels{10} = sprintf('sumdt');        if DATA.trigdt == 4 %template triggering - look at trig values            Labels{11} = sprintf('trigger');        else            Labels{11} = sprintf('sumdp',ispk);        end        Labels{12} = sprintf('sumdy',ispk);        Labels{14} = sprintf('MUsum',ispk);        Labels{15} = sprintf('sum2',ispk);        Labels{13} = sprintf('absdiff',ispk);        Labels{16} = sprintf('sum2-sum1');        Labels{17} = sprintf('sumP%d',DATA.probe(1)-1);        Labels{18} = sprintf('sumP%d',DATA.probe(1)+1);    end    