
mid=get_mid_number;

root=fullfile(topvfpaper,'PA_R21','VesselMask');
for i=6
    dname=sprintf('%s\\PVS%02d',root,i);
    if ~exist(dname,'dir')
      continue;
    end
     cd(dname);
    disp([i,j]);
    for j=1:2
        
        fphase=sprintf('vmask_Phase_PVS%02d_MID%d.mat',i,mid(i,j+1));
        
        if ~exist(fphase,'file')
            continue;
        end
      m1=load(fphase,'vmask');
      m2=load(sprintf('vmask_Mag_PVS%02d_MID%d.mat',i,mid(i,j+1)),'vmask');
      m1=clusterize2(m1.vmask);
      m2=clusterize2(m2.vmask);
      
      [m_ph,m_mag]=clusters_overlap(m1,m2);
      
      save(sprintf('vmask_Phase_omag_PVS%02d_MID%d',i,mid(i,j+1)),'m_ph');
      
    end

end

%% latest results saved in 7Jan2018, consistent with the paper.
%% data not ideal, do manual editting in IMG; data after editing saved in _22Jan2018.

