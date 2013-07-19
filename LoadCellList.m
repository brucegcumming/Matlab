function DATA = LoadCellList(DATA)
    cellfile = [DATA.name '/CellList.mat'];

    if exist(cellfile,'file')
        load(cellfile);
        nr = min([length(Expts) size(CellList,1)]);
        
%Dangerous to  change size of celllist. It will be written out again and so
%might overwrite datat that is wanted, just becuase and Expt is temporarily
%missing
%        DATA.CellList = CellList(1:nr,:);
        DATA.CellList = CellList;
        DATA.CellDetails = CellDetails;
        if ~isfield(CellDetails,'excludetrials')
            DATA.CellDetails.excludetrials = {};
        end
        if ~isfield(CellDetails,'exptids')
            for j = 1:size(Expts,1)      
                DATA.CellDetails.exptids(j) = Expts{j,1}.Header.exptno;
            end
        end
        if ~isfield(CellDetails,'trialids')
            for j = 1:size(Expts,1)
                DATA.CellDetails.trialids{j} = [Expts{j,1}.Trials.Trial];
            end
        else
            for j = 1:size(Expts,1)
                Tn = [Expts{j,1}.Trials.id];
                if length(setxor(DATA.CellDetails.trialids{j},Tn)) > 0 && j <= size(DATA.CellDetails.excludetrials,1)
                    fprintf('Expt Trial List Changed for Expt %d\n',Expts{j,1}.Header.exptno);
                    for k = 1:size(DATA.CellDetails.excludetrials,2)
                        if ~isempty(DATA.CellDetails.excludetrials{j,k,1})
                            [a, id] = ismember(Tn,DATA.CellDetails.excludetrials{j,k,1});
                        end
                    end
                end
            end
        end
        DATA = ConvertExclusion(DATA);
        if exist('CellListB','var')
            DATA.CellList(1:size(CellListB,1),1:size(CellListB,2),2) = CellListB;
        end
        if exist('CellChanges','var')
            DATA.CellChanges = CellChanges;
        else
            CellChanges = [];
        end
        if DATA.cellbackup == 0
        bakfile = strrep(cellfile,'.mat','back.mat');
        save(bakfile,'CellList','CellDetails','CellChanges');
        [a,b] = NetFilename(bakfile);
        if exist(b,'dir')
            fprintf('Backing up Cell List to %s\n',a);
            save(a,'CellList','CellDetails','CellChanges');
        else
            fprintf('Can''t Backup CellList to Network\n');
        end
        DATA.cellbackup = 1;
        set(DATA.toplevel,'UserData',DATA);
        end
    else
        DATA.CellList = zeros([length(DATA.exptid) DATA.nprobes 2])
    end


