function a=vec(a,roi,collapse_all)
if ~exist('roi','var')
    
a=a(:);

else
    
    if ~exist('collapse_all','var')
        collapse_all=true;
    end
    
    a=reshape(a,[length(roi(:)),length(a(:))/length(roi(:))]);
    
    a=a(roi(:)>0,:);
    
    if collapse_all
    a=a(:);
    end
end
