function m=mask_circle(dim,rad,center,include_equal)
% m=mask_circle(dim,rad,center,include_equal)
if ~exist('include_equal','var')
    include_equal=true;
end

if length(dim)==1
    dim=dim*ones(1,2);
end

[x,y]=meshgrid(1:dim(2),1:dim(1));

if include_equal
m=(x-center(2)).^2+(y-center(1)).^2<=rad^2;

else
    
m=(x-center(2)).^2+(y-center(1)).^2<rad^2;
end


