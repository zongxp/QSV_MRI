function fixGap(fname)

% assume that path{i} is 1:length(path{i})
a=load(fname);
i_ind_path=a.i_ind_path;

sz=size(a.c);
c=a.c;
ind=a.ind;
%% if there is a gap between the voxel and its path, delete that voxel

probclust=0*c;

nspath=0;
ri_ind_path=[];
for i=1:length(i_ind_path)
  if length(i_ind_path{i})==1
      ri_ind_path(end+1)=i;
      fprintf('Path %d has only 1 voxel\n',i);
      probclust(ind{i})=1;
      nspath=nspath+1;
  end
    
end

if nspath>0
i_ind_path(rpath)=[];
ind(rpath)=[];
fprintf('%d paths has only 1 voxel\n',nspath);
end
%%
probv=zeros(size(c));
for i=1:length(i_ind_path)
    
    if mod(i,50)==0
     %   disp(i);
    end
    pos=ind2subb(sz,ind{i});
    
    pospath=ind2subb(sz,ind{i}(i_ind_path{i}));
        
    while 1
        voxremoved=false;
      nvox=length(i_ind_path{i});
      while  nvox<length(ind{i})
            
            dist=sos(repmat(pos(nvox+1,:),[size(pospath,1),1])-pospath);
            
            
      
            if min(dist)>=2
                
                [tmp,indtmp]=min(dist);
                
                if hasgap(pos,nvox+1,i_ind_path{i}(indtmp))
                  
                  if ~voxremoved  
                    probclust(ind{i})=1;
                  end
                  
                  probv(ind{i}(nvox+1))=i;
                  c(ind{i}(nvox+1))=0;
                  ind{i}(nvox+1)=[];
                  
                  pos(nvox+1,:)=[]; 
                  voxremoved=true;
                  continue;
                end
            end
             nvox=nvox+1;
            
            
         end
           
         if ~voxremoved
             break;
         end
    end

end



if (sum(probv(:))>0)

    unqv=unique(probv(:));
  fprintf('%d PVSs found that have branches\n',  length(unqv)-1);
  
  fprintf('%d voxels will be removed to fix the problem\n',sum(probv(:)>0));
  %fprintf('%d voxels will be deleted from %d paths; paths:\n',sum(probv(:)>0),length(unique(probv(:)))-1);
   for i=1:length(unqv)
       if unqv(i)==0
           continue;
       end
      fprintf('path %d: %d voxels\n',unqv(i),sum(unqv(i)==probv(:))); 
   end
  
else
    disp('No braches found for any path.  Good!');
end

%disp(unique(probv(probv>0))');



c=zeros(sz);
pathroi=zeros(sz);

for i=1:length(i_ind_path)

    c(ind{i})=i;
    
    pathroi(ind{i}(i_ind_path{i}))=i;

end

fname=strtok(fname,'.');

voxsize=a.voxsize;

save([fname,'_Dv'],'ind','c','i_ind_path','probv','pathroi','probclust','voxsize');





    
function pos2=connect_pos(pos)

pos2=[];
for i=1:size(pos,1)-1
   
 n=pos(i+1,:)-pos(i,:);
 
 npix=ceil(sqrt(sum(n.^2)));
 
 for j=0:npix-1
    tmp=round(pos(i,:)+j*n/npix);
    if isempty(pos2) || any(tmp~=pos2(end,:))
      pos2(end+1,:) = tmp; 
    end
 end
 
    
end

if  any(pos(end,:)~=pos2(end,:))
      pos2(end+1,:) = pos(end,:); 
end


function res=hasgap(pos,i,j)

      pos2=connect_pos(pos([i,j],:));
            
      for i=2:size(pos2,1)-1
         
          if ~any(pos(:,1)==pos2(i,1) &pos(:,2)==pos2(i,2) & pos(:,3)==pos2(i,3))
            res= true;
            return;
          end
      end
      res=false;
      
      
      
      
      
      
                
            
function [maxlen,pathmax]=getCon2(pos, sz)

%%{
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
        
        
        
    

    
    



    
