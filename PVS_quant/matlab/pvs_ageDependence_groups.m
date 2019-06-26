function pvs_ageDependence_groups(roots,folders,ages)
% roots: 1*n
% folders: 1*n; each 1*x of directories; under which the *_Curv_Dv_stat.mat are
% ages: 1*n cell; each containing a vector of ages of the subjects in the
% group

sym={'mx','b^','cs','bo'};
    vol_thr=0.384;
    for j=1:length(roots)
        files_stat={};
        for i=1:length(folders{j})
            dir_str=dir(fullfile(roots{j},folders{j}{i},'*_Curv_Dv_stat.mat'));      
            files_stat{i}=fullfile(dir_str.folder,dir_str.name);        
        end
        pvs_age_dependence(files_stat,ages{j},vol_thr,sym{j});
    end
    
    