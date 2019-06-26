function d=Opac_mask(m,dim)
%m is a 3d mask


if dim==1
    v=[2,3,1];
elseif dim==2
    v=[1,3,2];
else   
    v=[1,2,3];
end

m=permute(m,v);

d=zeros(size(m(:,:,1)));

for i=1:size(m,1)
    for j=1:size(m,2)
        x=m(i,j,:);
        ind=find(x>0);
        if length(ind)>=1
          d(i,j)=m(i,j,ind(1));
        end
    end
end



