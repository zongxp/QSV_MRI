function res = diameter_vs_normalizedLength(diam,ind,matrix_size,vox_size,nseg)



pos=ind2subb(matrix_size,ind);

pos=pos.*repmat(vox_size,[length(ind),1]);

dist=sos(pos(1:end-1,:)-pos(2:end,:),2);

dist2=zeros(1,length(ind));
dist2(2:end)=dist;

dist2=dist2+mean(vox_size)/2;  %consider the end half voxel

dist2=cumsum(dist2);
dist_norm=dist2/(dist2(end)+mean(vox_size)/2);


gap=linspace(0,1,nseg);
    
    for j=1:size(gap,2)-1
        
        ind=dist_norm>=gap(j)&dist_norm<gap(j+1);
       
        res(j)=mean(diam(ind));
       
    end