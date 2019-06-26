function p=gui_pvs

p = parameter('PVS GUI');


p = add(p,'filename','prob files','Prob*.nii.gz'); 

p=add(p,'float','prob thr','0.8');
p=add(p,'float','clust size thr','6');
p = add(p,'button','gen mask','gui_pvs_callback(params)');


p = add(p,'filename','mask files','mask*.nii.gz'); 
p = add(p,'filename','T2 files','T2*.nii.gz'); 

p = add(p,'button','process','gui_pvs_callback(params)');


p = add(p,'button','gen bsub','gui_pvs_callback(params)');


%p=add(p,'float','clust thr (mm^3)','0.384');
% p = add(p,'directoryname','group 1','');
% 
% p = add(p,'directoryname','group 2','');
% 
% p = add(p,'directoryname','group 3','');

%p = add(p,'button','group analysis','edit group_ana.m');
% not so usefule; 
% p = add(p,'filename','T2 image','');
% p = add(p,'float','delay (s)',0.1);
% p = add(p,'float','max intensity',200);
% 
% p = add(p,'button','gen movie','gui_pvs_callback(params)');
p=add(p,'int','dimension','3');
p = add(p,'button','pvs mip','gui_pvs_callback(params)');
p = add(p,'button','pvs rotate gif','gui_pvs_callback(params)');

p = add(p,'button','Close','close');

if nargout==0
parametergui(p);
p=[];
end



