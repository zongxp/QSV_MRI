function [omask,nvox_sort]=clusterize2(mask,thr,ignore_diagnol)
%omask=clusterize2(mask,thr,ignore_diagnol)
% find connected clusters in the mask and remove voxels with cluster size
% less than thr.
% only version; NO clusterize or clusterize1.
%                    0: include all 
%if ignore_diagnol = 1: only ignore [1,1,1];
%                  = 2: ignore also [0,1,1], et al.,
%                  

if ~exist('ignore_diagnol','var')
    ignore_diagnol=0;
end

if ~exist('thr','var')
    thr = 1;
end

omask = zeros(size(mask));
sz = size(mask);

if length(sz)==2
    sz(3)=1;
end
    
clusters = {};
cluster_temp = zeros(sz(1)*sz(2)*sz(3),3);


ind=find(mask);
iclust=1;
for i=1:length(ind)
  ii=ind2subb(size(mask),ind(i));
  if length(ii)==2
      ii(3)=1;
  end
  if(mask(ii(1),ii(2),ii(3))==0) 
      continue;  % already in a cluster.
  end
  
  cluster_temp(1,:)=ii;
  
  
   mask(ii(1),ii(2),ii(3))=0;
         nvc=1; % number of voxels in a cluster.
            
         ivc = 1;
         while ivc<=nvc 
        
          nb = find_neighbors(cluster_temp(ivc,:),mask,ignore_diagnol);
        
          nnb = size(nb,1);   
         if nnb > 0 
          cluster_temp(nvc+1:nvc+nnb,:) = nb;
          nvc=nvc+nnb;
          for a=1:nnb
             mask(nb(a,1),nb(a,2),nb(a,3))=0;
          end
         end
        ivc = ivc+1;
         end
         
         
       clusters{end+1} = cluster_temp(1:nvc,:); 

       
end

nvox=zeros(1,length(clusters));

for i=1:length(clusters)
   nvox(i)=size(clusters{i},1); 
end

[nvox_sort,ind]=sort(nvox,'descend');

clusters=clusters(ind);

nclus = length(clusters);

for i=1:nclus
    if nvox_sort(i)<thr
        break;
    end
        
        for j=1:size(clusters{i},1)
                
            a=clusters{i}(j,:);
            omask(a(1),a(2),a(3)) =i;
        end
        
end
   %  fprintf('number of voxels in mask = %d\n',length(find(omask)));
    
        
        
function nb = find_neighbors(pos,mask,ignore_diagnol)


nb=[];
   sz = size(mask);
   if numel(sz) <3
       sz(3)=1;
   end
  
  for i=-1:1
      for j=-1:1
          for k=-1:1
              
              
            if(ignore_diagnol>0&&abs(i)+abs(j)+ abs(k)==3)
                continue;
            end
            
            if(ignore_diagnol>1&&abs(i)+abs(j)+abs(k)>=2)
                continue;
            end
            
            
              pos2=pos+[i,j,k];
              if pos2(1)<1 || pos2(1)>sz(1) ||pos2(2)<1 || pos2(2)>sz(2)||pos2(3)<1 || pos2(3)>sz(3)
                  continue;
              end
              
              if mask(pos2(1),pos2(2),pos2(3))>0
                  nb(end+1,:)=pos2;
              end
          end
      end
  end
      

    
   
     
   
