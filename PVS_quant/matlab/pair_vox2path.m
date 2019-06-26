function i_i_ind_path = pair_vox2path(ind,sz,i_ind_path)

% ind: indices for voxels in a vessel
% sz: image matrix size
% i_ind_path: indices for ind that defining the path.
% output:
% i_i_ind_path: should be the same length as the number of clusters.
% indices for i_ind_path that defines which voxel in the path each voxel
% belong to.

posc=ind2subb(sz,ind);

pos=ind2subb(sz,ind(i_ind_path));


for j=1:length(ind)
    
    dist=sos(repmat(posc(j,:),[size(pos,1),1])-pos);
    [mdist,itmp]=min(dist);
    i_i_ind_path(j)=itmp;
    
end
