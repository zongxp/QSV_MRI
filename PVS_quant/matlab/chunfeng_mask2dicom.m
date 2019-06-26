function m=chunfeng_mask2dicom(m,reverse)

%% convert data from mat file to correct orientation in nii format when created with make_nii, save_nii

if ~exist('reverse','var')
  reverse=false;
end

if   reverse
   m=flip(m,2);
   m=flip(m,1);
   m=permute(m,[2,1,3]);
   
else
m=permute(m,[2,1,3]);
%m=flip(m,3);
m=flip(m,1);

m=flip(m,2);
end