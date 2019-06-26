function gui_Slice2PVS


p = parameter('PVS GUI');

p = add(p,'filename','mask files','mask*.nii.gz'); 
p = add(p,'filename','T2 files','T2*.nii.gz'); 

p=add(p,'string','Slice center','');

p = add(p,'button','calculate','gui_pvs_callback(params)');


parametergui(p);
