function mask=niigz2mask_chunfeng(fname,per,flp)


gunzip(fname);

nii_name=strtok2(fname,'.');

a=load_nii(nii_name);

delete(nii_name);

b=a.img;
mask=uint8(b);
%a.img=flip(permute(mask,[2,1,3]),3);
%d=flip(flip(mask,1),2);
prefix=strtok(fname,'.');

if exist('per','var') && ~isempty(per)
    mask=permute(mask,per);
end

if exist('flp','var')
   for i=1:length(flp)
      mask=flip(mask,flp(i)); 
   end
    
end

if nargout==0
  save(prefix,'mask');
end


