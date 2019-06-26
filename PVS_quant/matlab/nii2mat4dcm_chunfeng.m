function nii2mat4dcm_chunfeng(fname)


a=load_nii(fname);

d=flip(permute(a.img,[2,1,3]),3);

prefix=strtok(fname,'.');

save([prefix,'.mat'],'d');

