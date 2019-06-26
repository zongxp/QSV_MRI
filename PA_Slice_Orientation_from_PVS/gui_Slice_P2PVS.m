function gui_Slice_P2PVS(dname,simple)

if ~exist('simple','var')
    simple=true;
end

if ~exist('dname','var')
    dname='';
end

p = parameter('PVS GUI');

if simple
p = add(p,'directoryname','dicom dir name',dname);

p=add(p,'string','Slice center','');

p = add(p,'button','calculate','gui_Slice_P2PVS_callback(params,1)');
else    
p = add(p,'filename','mask files','mask*.nii.gz');
p = add(p,'filename','T2 files','T2*.nii.gz');
p = add(p,'filename','orient files','mask*.nii.gz');

p=add(p,'string','Slice center','');

p = add(p,'button','calculate','gui_Slice2PVS_callback(params,0)');

end
parametergui(p);