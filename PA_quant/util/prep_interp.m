function [v_pc,s_pc,s_static]=prep_interp(TR, T1, T2, TE, thk, FA)

f_prep=['TR_',num2str(TR),'_T1_',num2str(T1),'_T2_',num2str(T2),'_TE_',num2str(TE),'_thk_',num2str(thk),'_FA_',num2str(FA),'.mat'];
f_prep=strrep(f_prep,' ','_');
for i=1:10
    f_prep=strrep(f_prep,'__','_');
end


root=fileparts(mfilename('fullpath'));
dname=fullfile(root,'prep_interp');

if ~exist(dname,'dir')
    mkdir(dname);
end

d_prep=fullfile(dname,f_prep);

if ~exist(d_prep,'file')
    [v_pc,s_pc,s_static]=calc_FlowEnhancement(TR,T1,TE,T2,thk,FA,'sinc');
    save(d_prep,'v_pc','s_pc','s_static');
else
    load(d_prep);
end

