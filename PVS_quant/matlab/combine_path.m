function combine_path(fname)


a=load(fname);
path=a.path;

sz=size(a.c);

c=a.c;
ind=a.ind;
%% connect paths if they are close to each other 
rpath=[];
for i=1:length(path)
  if length(path{i})==1
      rpath(end+1)=i;
      fprintf('Path %d has only 1 voxel\n',i);
      
  end
    
end


path(rpath)=[];
ind(rpath)=[];


pcombine=[];
for i=1:length(path)
  
      
    if mod(i,50)==0
        disp(i);
    end
    
    pos1=ind2subb(sz,ind{i}(path{i}));
  
    for j=2:length(path)
     
     if j<=i
            continue;
     end
   
    
        
       pos2=ind2subb(sz,ind{j}(path{j}));
    
        if sos(pos1(1,:)-pos2(1,:))<3
            
            dpos1=pos1(1,:)-pos1(2,:);
            dpos2=pos2(2,:)-pos2(1,:);
            
        elseif   sos(pos1(1,:)-pos2(end,:))<3
            
            dpos1=pos1(1,:)-pos1(2,:);
            dpos2=pos2(end-1,:)-pos2(end,:);
        elseif sos(pos1(end,:)-pos2(1,:))<3
            
            dpos1=pos1(end,:)-pos1(end-1,:);
            dpos2=pos2(2,:)-pos2(1,:);
            
        elseif sos(pos1(end,:)-pos2(end,:))<3
            
            dpos1=pos1(end,:)-pos1(end-1,:);
            dpos2=pos2(end-1,:)-pos2(end,:);
            
        else
            continue;
        end
        
         if sum(dpos1.*dpos2)/sos(dpos1)/sos(dpos2)>sqrt(3)/2
                  pcombine(end+1,:)=[i,j];             
          end
        
        
    end
    
    
end


    fprintf('%d pathes will be combined; paths:\n',length(unique(pcombine(:))));
    disp(pcombine);
    upcombine=unique(pcombine(:));
  if length(upcombine)<length(pcombine(:))
    for i=1:length(unique(pcombine(:)))
      
        if length(upcombine(i)==pcombine(:))==1
            continue;
        end
        
        for j=1:size(pcombine,1)
            
           if pcombine(j,2)==upcombine(i)
               
               pcombine(j,:)=fliplr(pcombine(j,:));
           end
           
            
        end
        
        
    end
  end

rpathc=[];
save temp

load temp;
newpath=zeros(size(c));
for i=1:size(pcombine,1)
    
    icl2=pcombine(i,2);
    
    icl1=pcombine(i,1);
   c(c==icl2)=icl1;
   
    pos1=ind2subb(sz,ind{icl1}(path{icl1}));
    pos2=ind2subb(sz,ind{icl2}(path{icl2}));
    
        if sos(pos1(1,:)-pos2(1,:))<4
            dlen= sos(pos1(1,:)-pos2(1,:));
            posnew=connect_pos([pos1(1,:);pos2(1,:)]);
            
        elseif   sos(pos1(1,:)-pos2(end,:))<4
             dlen= sos(pos1(1,:)-pos2(end,:));
            posnew=connect_pos([pos1(1,:);pos2(end,:)]);
            
        elseif sos(pos1(end,:)-pos2(1,:))<4
            
            dlen= sos(pos1(end,:)-pos2(1,:));
            posnew=connect_pos([pos1(end,:);pos2(1,:)]);
            
        else
            
            dlen= sos(pos1(end,:)-pos2(end,:));
            posnew=connect_pos([pos1(end,:);pos2(end,:)]);
        end
   
   
   for j=1:size(posnew,1)
       
       c(posnew(j,1),posnew(j,2),posnew(j,3))=icl1;
   end
   
   
    rpathc(end+1)=icl2;
    
    ind{icl1}=find(c(:)==icl1);   
    

    [path{icl1},ind{icl1}]=findPathMatlab(ind{icl1},size(c));
     newpath(ind{icl1}(path{icl1}))=icl1;
end



path(rpathc)=[];
ind(rpathc)=[];


c=zeros(sz);
pathroi=zeros(sz);

for i=1:length(path)   
    pathroi(ind{i}(path{i}))=i;
    c(ind{i})=i;
end

prefix=strtok(fname,'.');



save([prefix,'_Comb'],'pathroi','c','ind','path','newpath');





    
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
      
      
      
      
      
      
function [pathmax,ind,maxlen]=findPathMatlab(ind, sz)

%% pos is the ind of voxels in the cluster.1*nvox
%% sz: matrix size.
%%{


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
    
    
    
    
    
    
    
    
        
        
    
    
    



    
