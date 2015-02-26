function DATA = SetCluster(DATA, eid, probes, varargin)


if isfield(DATA,'AllClusters')
    return;
end
id = ismember(probes,1:size(DATA.cluster,2));
probes = probes(id);
for j = probes(:)'
        pdef(j) = 0;
        ndef(j) = 0;
%copy any clusters defined ine Expts{eid} to DATA.cluster        
        for k = 1:size(DATA.cluster,1)
            if isfield(DATA.Expts{eid},'Cluster') && cmb.iscluster(DATA.Expts{eid}.Cluster,k,j)
                DATA.cluster{k,j} = DATA.Expts{DATA.currentexpt(1)}.Cluster{1,j};
                DATA.cluster{k,j}.touched = 1;
                ndef(j) = ndef(j)+1;
                pdef(j) = pdef(j)+1;
            else
                DATA.cluster{k,j}.touched = 0;
            end
        end
        
        if ndef(j) == 0 
%get number of expt associated with current cluster            
            if isfield(DATA.cluster{k,j},'exptno')
                cleid = DATA.cluster{k,j}.exptno;
            else
                cleid = 0;
            end
            if cmb.iscluster(DATA.cluster,1,j) & DATA.state.applylastcluster
                fprintf('Applying Cluster From Expt %d to probe %d\n',cleid,j)
                for k = 1:size(DATA.cluster,1)
                    if cmb.iscluster(DATA.cluster,k,j)
                        DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j} = DATA.cluster{k,j};
                            DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j}.fromlast = cleid;
                        ndef(j) = ndef(j)+1;
                    end
                end
            elseif DATA.plot.autoclustermode < 3
                DATA = cmb.AutoCut(DATA, eid);
            end
            if ~isfield(DATA.cluster{1,j},'exptno')
                DATA.cluster{1,j}.exptno = 1;
            end
        end
    end
