function run_vmask_medianFilter_phase(alpha_thr)
mid=get_mid_number;
q=get_quality;

if ~exist('alpha_thr','var')
    alpha_thr=0.05;
end

for i=6
    for j=1:2
     
        savename=sprintf('PVS%02d/vmask_Phase_PVS%02d_MID%d.mat',i,i,mid(i,j+1));
        if exist(savename,'file')
            continue;
        end   
        
        fname=sprintf('../mat_phase_mag/PVS%02d_mat/Phase_PVS%02d_MID%d.mat',i,i,mid(i,j+1));
        if ~exist(fname,'file')
         continue;
        end
        
     disp([i,j]);
         pc=ri_d1(fname,'','','pc');
        
    params.data=double(pc)*pi/18000;
    params.patchSize=round(10/0.3125*2);
    
    fname2=name4pat(sprintf('../MagSNRMap/PVS%02d/SNRMap_Mag_MID%d*.mat',i,mid(i,j+1)));
    params.SNR_mag=ri(fname2);
 
    iscan=get_slice_index_T1(i,j);
    fname=name4pat(sprintf('../wm_mask/wm_mask_shrink/wm_mask_scan%d_PVS%02d_shk_nohole.mat',iscan,i));
    
    params.ispc=true;
    params.wm_mask=ri_d1(fname);
    params.interp_factor=[2,3.3333];
    params.alpha_2tail=alpha_thr;
    vmask=vmask_medianFilter(params);
    
    % params contain:
    % pc: map of PC; assume positive value for vessels; in radians;
    % patchSize: size for meidan filter;
    % SNR_mag: magnitude SNR map;
    % mask: white matter mask
    % interp_factor: interpolation factor needed for Bonferroni correction.
    % alpha_2tail: default 0.05;
    if ~exist(sprintf('PVS%02d',i),'dir')
        mkdir(sprintf('PVS%02d',i));
    end
     
    save(savename,'vmask');
    end
end

%% latest results saved in 7Jan2018, consistent with the paper.
%% data not ideal, do manual editting in IMG; data after editing saved in _22Jan2018.


function res=get_slice_index_T1(isub,iscan)

loc=get_slice_locations;

res=find(loc(isub,iscan+1)==loc(isub,2:iscan+1),1,'first');

