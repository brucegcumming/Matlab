function res = RunExptCmds(name,expts, varargin)
%res =  RunExptCmds(name,expts, varargin)
%simple version of RunnAllVPcs. Builds FullV name for each number in 
%expts. Runs AllVPcs(expt, varargin{:}).

if iscell(name)
    CheckErrors(name);
    return;
end
for j = 1:length(expts)
    filename = [name '/Expt' num2str(expts(j)) 'FullV.mat']
    res.cls{j} = AllVPcs(filename, varargin{:});
end



function CheckErrors(X)

for j = 1:length(X)
    for k = 1:length(X{j})
        C = X{j}{k};
        for e = 1:length(C.errs)
            fprintf('%s\n',C.errs{e});
        end
        for e = 1:length(C.errstates)
            fprintf('%s\n',C.errstates{e}.message);
        end
    end
end

