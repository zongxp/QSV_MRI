function gui_Slice_P2PVS_callback(p,simple)


if simple
   dname = get(p,'dicom dir name');
   fc=fullfile(sprintf('mask_%s',dname),sprintf('mask_%s_PVS.mat',dname));
   fo=fullfile(sprintf('mask_%s',dname),sprintf('mask_%s_PVS_orient.mat',dname));
   f_T2=sprintf('%s.nii.gz',dname);
  
   center=get(p,'Slice center');
   center=get_center(center);
   prefix=sprintf('Slice_%s_%d_%d_%d',dname,round(center));
   
else
    
    
end


Slice_P2PVS(fo,fc,prefix,center,f_T2);


function res=get_center(center)

res=str2num(center);
if any(isnan(res))
    
    
end