function prefix=gen_pvsmask(prob_file,prob_thr,clust_size)

if ~exist('prob_thr','var')
    prob_thr=0.8;
end

if ~exist('clust_size','var')
    clust_size=6;
end

nii_name=strtok2(prob_file,'.');
prefix=strtok2(nii_name,'.');
prefix=filename_prefix(prefix(6:end),'mask_');

if ~exist([prefix,'.nii.gz'],'file')
    d=load_untouch_niigz(prob_file);
    
    [m,nvox_sort]=clusterize2(d.img>prob_thr,clust_size);
    
    d.img=m>0;
    fprintf('%s - Number of PVSs found = %d\n',filename(prefix),max(m(:)));
    save_untouch_niigz(d,prefix,true);
end