function dcm2pvs(dname)
th=tic;
pvsseg_remote(dname);


prefix=gen_pvsmask(sprintf('prob_%s.nii.gz',dname));
process_pvs([prefix,'.nii.gz'],[],1);

gui_Slice_P2PVS(dname);
toc(th);

% %%
% h=img;
% par_fname.ufile=[dname,'.nii.gz'];
% par_fname.uroifile=fullfile(prefix,[prefix,'_PVS.mat']);
% 
% par_vname.uroifile='c';
% 
% 
% set_data_img(par_fname,par_vname,h);
% a=findobj(h,'Tag','first_slice');
% set(a,'string','103');
% 
% showfull(h);
% swapxz(h);
% swapxy(h);
% fliplr(h);
% flipud(h);


%%



