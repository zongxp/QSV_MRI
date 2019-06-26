function y=repmat2(c_image,nd)
% y=repmat2(c_image,nd); replicate each voxel nd times 
% directions 
% nd can be 1, 2, or 3 elements;
% when has only 1 element, equivelant to [nd,nd,1];
% when has only 2 elements, equivelant to [nd,1];

if length(nd)==1
    nd=[nd,nd,1];
elseif length(nd)==2
    nd=[nd,1];
end

y=zeros(nd(1),size(c_image,1),nd(2),size(c_image,2),nd(3),size(c_image,3),size(c_image,4));

for l=1:size(c_image,4)
for k=1:size(c_image,3)
    for i=1:size(c_image,1)
        for j=1:size(c_image,2)
            ind=c_image(i,j,k,l);
            
         y(:,i,:,j,:,k,l)=ind;
        end
    end
    
end
end

y=reshape(y,[size(c_image,1)*nd(1),size(c_image,2)*nd(2),size(c_image,3)*nd(3),size(c_image,4)]);
    
    
    
    
    