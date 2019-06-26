function nii2mask_chunfeng(fname)


a=load_nii(fname);

mask=uint8(a.img>0.5);
%a.img=flip(permute(mask,[2,1,3]),3);
a.img=flip(flip(mask,1),2);  %why do I need this?

prefix=strtok(fname,'.');

[id,suf]=strtok2(prefix,'_');

a.hdr.dime.bitpix = 8;
a.hdr.dime.datatype=2;

save_nii(a,['mask_',suf(2:end),'.nii']);

