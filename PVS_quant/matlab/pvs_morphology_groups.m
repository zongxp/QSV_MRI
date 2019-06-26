function pvs_morphology_groups(roots,folders,labels)
%  pvs_morphology_groups(roots,folders,labels)
% roots: 1*n cell
% folders: 1*n cell; each 1*x of directories; under which the *_Curv_Dv_stat.mat are
% 1*n cell

sym={'mx-','c^-','ks-','bo-'};
    vol_thr=0.384;
    for j=1:length(roots)
        files_stat={};
        
        for i=1:length(folders{j})
       
            dir_str=dir(fullfile(roots{j},folders{j}{i},'*_Curv_Dv_stat.mat'));      
            files_stat{i}=fullfile(dir_str.folder,dir_str.name);        
        end
        pvs_stat(files_stat, vol_thr,sym{j});
    end
    
    subplot(2,2,1);
    legend(labels);