function [pathmax,ind,maxlen]=findPath(ind, sz)

%% pos is the ind of voxels in the cluster.1*nvox
%% sz: matrix size.
%%{



pos=ind;

m=zeros(sz);
m(ind)=1;

c=clusterize2(m,2);

if max(c(:))>1
    disp('Any not find a path for more than two clusters');
    return;
end
[len,pair,path,nfound]=find_neighbors2(pos,sz);


n=length(pos);
tic;
while sum(nfound)<n*(n-1)/2
   disp(sum(nfound));
  for i=1:n
    
    for j=2:n
        if (i>=j)
            continue;
        end
        
        if any(j==pair(i,:)) 
            continue;
        end
            
        [lentmp,path_tmp]=share_common(n,i,j,pair,len,path,nfound);
        
         if ~isinf(lentmp)
             
             
             nfound(i)=nfound(i)+1;
             pair(i,nfound(i))=j;
             len(i,nfound(i))=lentmp;   
             path{i,nfound(i)}=path_tmp;
         end      
    end
  end
  
end

  maxlen=0;
  for i=1:n
   [lentmp,ind]=max(len(i,:));
   if maxlen<lentmp
    pathmax=path{i,ind};
    maxlen=lentmp;
   end
  end
  
   xpath=setdiff(1:length(ind),pathmax);
     
       ind=ind(:);
    ind=[ind(pathmax);ind(xpath)];
 
    pathmax=1:length(pathmax);

  
  
    function [minlen,path_new]=share_common(n,i,j,pair,len,path,nfound)
        
        minlen=Inf;
        path_new=[];
     
    %{   
       if isempty(intersect(pair(i,:),pair(j,:)))
           return;
       end
     %}
     
        for k=1:n
            
            if k==i || k==j 
                continue;
            end
            
            i1=i;
            k1=k;
            if k<i
             i1=k;
             k1=i;
            end
            
            ii1= find(pair(i1,1:nfound(i1))==k1);
            
            
            if isempty(ii1)
                continue;
            end
            
            
            i2=j;
            k2=k;
            if k<j
             i2=k;
             k2=j;
            end
            
           
            
            ii2= find(pair(i2,1:nfound(i2))==k2);
           
            if isempty(ii2)
                continue;
            end
               if minlen>len(i1,ii1)+len(i2,ii2)
                  
                       p1=path{i1,ii1};
                       p2=path{i2,ii2};
                 
                  
                     if p1(end)==p2(end)  
                      path_new=[p1,fliplr(p2(1:end-1))];
                     elseif p1(end)==p2(1)
                      path_new=[p1,p2(2:end)];
                     elseif p2(end)==p1(1)  
                      path_new=[p2,p1(2:end)];
                     elseif p2(1)==p1(1)
                         path_new=[fliplr(p1),p2(2:end)]; 
                     else
                         error('this should not happen');
                     end
                  %}   
                       minlen=len(i1,ii1)+len(i2,ii2);
                  
                        
               end
               
               
            
            
            
        end
        
        
        
function [len,pair,path] = find_neighbors(pos,sz)

    len=[];
    pair=[];
    path={};
    for i=1:length(pos)
       pos1=ind2subb(sz,pos(i));
    
        for j=2:length(pos)
            
          if i>=j
              continue;
          end
          
          pos2=ind2subb(sz,pos(j));
          
          if ~any(abs(pos1-pos2)>1)
              
              pair(end+1)=pairIndex(i,j);
              len(end+1)= sqrt(sum(abs(pos1-pos2).^2));
              
              path{end+1}=[i,j];
          end
  
        end
    end
    
    function [len,pair,path,nfound] = find_neighbors2(pos,sz)

        n=length(pos);
        
    len=zeros(n,n);
    pair=zeros(n,n);
    nfound=zeros(1,n);
    path=cell(n,n);
    
    for i=1:length(pos)
       pos1=ind2subb(sz,pos(i));
    
        for j=2:length(pos)
            
          if i>=j
              continue;
          end
          
          pos2=ind2subb(sz,pos(j));
          
          if ~any(abs(pos1-pos2)>1)
          
          nfound(i)=nfound(i)+1;    
          pair(i,nfound(i))=j;
          
              len(i,nfound(i))= sqrt(sum(abs(pos1-pos2).^2));
              
              path{i,nfound(i)}=[i,j];
          end
  
        end
    end
    
    function ind=pairIndex(i,j)
        
        if (i>j)
            tmp=j;
            j=i;
            i=tmp;
        end
        ind=i*10000+j;
        
        
        
    