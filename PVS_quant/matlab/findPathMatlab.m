function [pathmax,ind,maxlen]=findPathMatlab(ind, sz)

%% pos is the ind of voxels in the cluster.1*nvox
%% sz: matrix size.
%% use mathlab graphshortestpath function to find the path 


pos=ind2subb(sz,ind);

for i=1:3
[tmp,imax(i)]=max(pos(:,i));
[tmp,imin(i)]=min(pos(:,i));
end

iall=[imax,imin];
dist=0;
pair=[];
for i=1:6
    for j=2:6
        
        if j<=i 
            continue;
        end
        
        if iall(i)==iall(j)
            continue;
        end
        
        disttmp=sos(pos(iall(i),:)-pos(iall(j),:));
        if dist<disttmp
            pair=[iall(i),iall(j)];
            dist=disttmp;
        end
        
    end
end


[s,w]=find_neighbors2(pos);

sp=sparse(s(:,1)',s(:,2)',w);
sp(size(sp,1)+1:size(sp,2),:)=0;

sp=tril(sp + sp');
[maxlen,pathmax] = graphshortestpath(sp,pair(1),pair(2),'Directed',false);

ind=ind(:);


xpath=setdiff(1:length(ind),pathmax);
ind=[ind(pathmax);ind(xpath)];
pathmax=1:length(pathmax);




        
        

    
    function [p,w] = find_neighbors2(pos)

       
    p=[];
    w=[];
    for i=1:size(pos,1)
       pos1=pos(i,:);
    
        for j=2:size(pos,1)
            
          if i>=j
              continue;
          end
          
          pos2=pos(j,:);
          
          if ~any(abs(pos1-pos2)>1)
          
             p(end+1,:)=[i,j];
             
                
             w(end+1)= sqrt(sum(abs(pos1-pos2).^2));
              
             
          end
  
        end
    end
    
    
    
    
    
    
    
    
        
        
    