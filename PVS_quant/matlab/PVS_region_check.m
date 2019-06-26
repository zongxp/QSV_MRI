function [ind_n,ind_n2,roi_label]=PVS_region_check(pvs_seg_fname, roi,suffix)
% pvs_seg_fname: output from pvs_quantify
% roi: 4d matrix containing size(roi,4) rois

%ind_n: PVS indices
%ind_n2: PVS indices excluding ambiguous PVS


if ~exist('suffix','var')
    suffix='_regSpc';
end
[ind_n,ind_n2,roi_label,roi_nvox] = get_regions(pvs_seg_fname,roi);


l=ri(pvs_seg_fname,'','','l');
v=ri(pvs_seg_fname,'','','v');
dpath=ri(pvs_seg_fname,'','','dpath');
vsize=ri(pvs_seg_fname,'','','voxsize');
fname=filename_append(pvs_seg_fname,suffix);


save(fname, 'roi','ind_n2','ind_n', 'l', 'v', 'dpath','roi_label','roi_nvox');
copy_orient(pvs_seg_fname,fname);

% else
%   fname=filename_append(pvs_seg_fname,suffix);
%     
%  mat_addv(fname,'roi_nvox',roi_nvox);
% end

function res=calc_roi_nvox(roi) 

for i=1:size(roi,4)
   res(i)=sum(vec(roi(:,:,:,i)>0));
end


function  [ind_n,ind_n2,roi_label,roi_nvox] = get_regions(pvs_seg_fname,roi)



roi_nvox=calc_roi_nvox(roi);

a=load(pvs_seg_fname);
path=a.path;
ind=a.ind;

nroi=size(roi,4);
ind_n=cell(1,nroi);  % PVS indices 
ind_n2=cell(1,nroi);  % PVS indices excluding ambiguous PVS
roi_label = 0*roi(:,:,:,1);
sz=size(roi);
roi=reshape(roi,[prod(sz(1:3)),sz(4)]);


for i=1:length(path)
 
    indp=ind{i}(path{i});
    
    ntmp=sum(roi(indp,:),1);
    
    if ~any(ntmp>0)
        continue;
    elseif sum(ntmp>0)==1  % PVS only intersect one region
        ind_n{ntmp>0}(end+1)=i;
        ind_n2{ntmp>0}(end+1)=i;
        roi_label(ind{i})=find(ntmp>0);
    else % PVS intersect more than one region
        
        [nmax,ireg]=max(ntmp);  
        ind_n{ireg}(end+1)=i;
        roi_label(ind{i})=ireg;
    end  
    
end

disp(' ');
