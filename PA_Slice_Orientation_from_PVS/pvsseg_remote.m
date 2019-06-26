function pvsseg_remote(dname)

dcm2nii(dname);

a=load(sprintf('DCMPar_%s.mat',dname));
if a.par.VoxelSize(1)<0.39  %interpolated
    undersample_image(sprintf('%s.nii.gz',dname),[2,2,1]);
    movefile(sprintf('%s_ds.nii.gz',dname),sprintf('%s.nii.gz',dname));
end


upload_data(dname);

run_pvsseg(dname);

download_prob(dname);



function upload_data(dname)

cmd='winscp.exe /command "open erwin" "option batch on" "cd /data/zong/PVS"';
cmd=sprintf('%s "lcd "%s""',cmd,pwd);
cmd=sprintf('%s "put %s.nii.gz"',cmd,dname);
fail_free(cmd);




function run_pvsseg(dname)

cmd='winscp.exe /command "open erwin" "option batch on" "cd /data/zong/PVS"';
cmd=sprintf('%s "lcd "%s""',cmd,pwd);
cmd=sprintf('%s "call pvsseg.py %s.nii.gz"',cmd,dname);
fail_free(cmd);


function download_prob(dname)

cmd='winscp.exe /command "open erwin" "option batch on" "cd /data/zong/PVS"';
cmd=sprintf('%s "lcd "%s""',cmd,pwd);
cmd=sprintf('%s "get prob_%s.nii.gz"',cmd,dname);
fail_free(cmd);
