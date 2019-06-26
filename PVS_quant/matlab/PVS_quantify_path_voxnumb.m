function PVS_quantify_path_voxnumb
%fname = 'PVSLength_testzb_Curv_Dv_Comb_Dv.mat';
%fname='PVSLength_HVB238b_Curv_Dv_Comb_Dv.mat';
%fname='PVSLength_HVB239c_Curv_Dv_Comb_Dv.mat';

d=236:241;

for i=1:length(d)
   
fname=sprintf('PVSLength_FinalLabel_%d_Curv_Dv_stat',d(i));

a=load(fname);
path=a.path;

nvox=zeros(1,length(path));

for j=1:length(path)
 
nvox(j) = length(path{j});        
    
    
end

disp(max(nvox));

end