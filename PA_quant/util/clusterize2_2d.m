function mask=clusterize2_2d(mask,thr,ignore_diagnol)
%[mask,clusters]=clusterize2(mask,thr)
% find connected clusters in the mask and remove voxels with cluster size
% less than thr.
% only version; NO clusterize or clusterize1.

if ~exist('ignore_diagnol','var')

    ignore_diagnol=false;
end
if ~exist('thr','var')
    thr = 1;
end

for i=1:size(mask,3)
    
    tmp=mask(:,:,i);
    if sum(tmp(:)>0)<thr
        mask(:,:,i)=0;
        continue;
    end
    
    mask(:,:,i)=clusterize2(mask(:,:,i),thr,ignore_diagnol);
    
    
    if i>1
        tmp=mask(:,:,i);
        tmp(tmp>0)=tmp(tmp>0)+max(vec(mask(:,:,1:i-1)));        
       %mask(:,:,i)=mask(:,:,i)+max(vec(mask(:,:,1:i-1)));
        mask(:,:,i)=tmp;
    end
    
end


