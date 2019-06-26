function [val,eval,val_voxel]=mean_roi(data,mask)
% [val,nv]=mean_roi(data,mask)
% do average over roi defined in mask.
% data can be file name or 4d matrix
% mask can be file name or 3d matrix

    if isa(data,'char')
        data = ri(data);
    end
    
    if isa(mask,'char')
        mask = ri(mask);
    end

    data=double(data);
    nd=ndims(data);
    nm=ndims(mask);
    lm=length(mask(:));
    ld=length(data(:));
    
    data2=reshape(data,[lm,ld/lm]);
    
    val_voxel=data2(mask(:)>0,:);
    
    data3=mean(val_voxel,1);
    edata3=std(val_voxel,[],1);
    
    sz=size(data);
    rshpsz=sz(nm+1:end);
    if length(rshpsz)>1
     val=reshape(data3,rshpsz);
     eval=reshape(edata3,rshpsz);
    else
        val=data3(:);
        eval=data3(:);
    end
    
    nv = length(find(mask>0));
    
%disp([mfilename ' finish in ', num2str(toc), ' s']);