function [m,pos2]=Mask_between_2voxels(pos,sz)
%pos: 2*3 the [x,y,z] coordinates of two voxels

m=zeros(sz);

 n=pos(2,:)-pos(1,:);
 
 npix=ceil(sqrt(sum(n.^2)));
 pos2=[];
 
 for j=0:npix
    tmp=round(pos(1,:)+j*n/npix);
      
      if isempty(pos2)|| any(pos2(end,:)~=tmp)
          pos2(end+1,:)=tmp;
          m(tmp(1),tmp(2),tmp(3))=1; 
 
      end
          
 end
 
