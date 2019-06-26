function checkPaths(fname,deletepath)


a=load(fname);
path=a.path;

sz=size(a.c);
l=a.l;
v=a.v;
r=a.r;
c=a.c;
%% 

for i=1:length(path)
   ind{i}=find(c(:)==i);
   indb{i}=ind{i};  % new ind after changes to the voxels
   if mod(i,50)==0
       disp(i);
   end
end
%{
pathroi=zeros(size(c));
for i=1:length(path)
    
    pathroi(ind{i}(path{i}))=i;
end
roi=pathroi;
fname=strtok(fname,'.');
save([fname,'_pathroi'],'roi');
%}
%% if there is a gap between the voxel and its path, delete that voxel

probv=zeros(size(c));
for i=1:length(path)
    
    if mod(i,50)==0
        disp(i);
    end
    pos=ind2subb(sz,ind{i});
    
    pospath=ind2subb(sz,ind{i}(path{i}));
        
    doclust=false;
         for k1=1:size(pos,1)
            
            dist=sos(repmat(pos(k1,:),[size(pospath,1),1])-pospath);
            
            
      
            if min(dist)>=2
                
                [tmp,indtmp]=min(dist);
                
                if hasgap(pos,k1,path{i}(indtmp))
                  probv(ind{i}(k1))=i;
                  c(ind{i}(k1))=0;
                  doclust=true;
                end
                
            end
             
         end
           
         if doclust          
          tmproi=clusterize2(c==i);
          
          c(tmproi>1)=0;
          probv(tmproi>1)=i;
          indb{i}=find(c(:)==i);
         end
         
    

end

fprintf('%d voxels will be deleted from %d paths; paths:\n',sum(probv(:)>0),length(unique(probv(:)))-1);
disp(unique(probv(probv>0)));

roi=probv;
fname=strtok(fname,'.');
save([fname,'_deletevox']);

%% if two paths are too close, remove the smaller one.
%%{
%load all

%{
path(rpath)=[];
l(rpath)=[];
v(rpath)=[];
r(rpath)=[];
ind(rpath)=[];
fname=strtok(fname,'.');
roi=clusterize2(deletec);
save([fname,'_deleted_path'],'roi');


%}
%% connect paths if they are close to each other 


pcombine=[];
for i=1:length(path)

    
    if mod(i,50)==0
        disp(i);
    end
    
    pos1=ind2subb(sz,ind{i}(path{i}));
  %{      
    if any(i==rpath)  
        continue;
    end
    %}
    for j=2:length(path)
     
     if j<=i
            continue;
     end
    %{ 
    if any(j==rpath) 
        continue;
    end
    %}
    
        
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
            
        elseif sos(pos1(end,:)-pos2(end,:))<=3
            
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
newpath=zeros(size(c));
rpathc=[];
save temp
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
   
   l(icl1)=l(icl1)+l(icl2)+dlen;
   
   for j=1:size(posnew,1)
       
       c(posnew(j,1),posnew(j,2),posnew(j,3))=icl1;
   end
   v(icl1)=sum(c(:)==icl1);
   r(icl1)=sqrt(v(icl1)/l(icl1)/pi);
    rpathc(end+1)=icl2;
    
  %  tmproi=clusterize2(newpath==icl1,2);
  %  ind{icl1}=find(tmproi==1);
  
    ind{icl1}=find(c(:)==icl1);   
    indb{icl1}=ind{icl1};
    newpath(ind{icl1})=icl1;
  
    [tmp,path{icl1}]=getCon2(ind{icl1},size(c));
    
end

prefix=strtok(fname,'.');
roi=newpath;
save([prefix,'_combined_path']);

%% delete some paths if too close to another path
if ~deletepath
    roi=clusterize2(c);
    save([prefix,'_corrected_nodelpath']);
    
    return;
end

    
rpath=[];
%%{
deletec=zeros(size(c));
for i=1:length(path)

    
    if mod(i,50)==0
        disp(i);
    end
     pos1=ind2subb(sz,ind{i}(path{i}));
     %pos1=ind2subb(sz,ind{i});
        
    if any(i==rpath) || any(i==rpathc)
        continue;
    end
    for j=2:length(path)
     
         
    if any(j==rpath) || any(j==rpathc)
        continue;
    end
    
     if j<=i
            continue;
     end
        
     %  pos2=ind2subb(sz,ind{j}(path{j}));
      pos2=ind2subb(sz,indb{j});
    
       
        %%{
    
        if sos(mean(pos1,1)-mean(pos2,1))>l(i)+l(j)
            continue;
        end
        
        
        lmin=i;
        if (l(j)<l(i))
            lmin=j;
            
            pos1=ind2subb(sz,ind{j}(path{j}));
            pos2=ind2subb(sz,indb{i});
            
            
        end
        nn=0;
        
           for k1=1:size(pos1,1)
            
            dist=sos(repmat(pos1(k1,:),[size(pos2,1),1])-pos2);
            
            if min(dist)<3
                nn=nn+1;

            end
            
           end
      
        
        if  nn/l(lmin)>=0.3
            
          
            rpath(end+1)=lmin;
            
            
            deletec(c==lmin)=lmin;
            c(c==lmin)=0;
        end
        
        
    end
    
    
end

fprintf('%d paths will be removed; path: \n',length(rpath));
disp(rpath);
%}

path([rpath,rpathc])=[];
v([rpath,rpathc])=[];
r([rpath,rpathc])=[];
l([rpath,rpathc])=[];
roi=clusterize2(c);
ind([rpath,rpathc])=[];
save([prefix,'_corrected_delpath']);
%{
figure;subplot(1,3,1);
hist(l,0:2:160);

set(gca,'FontSize',12);
xlabel('PVS length (Unit: Voxel)');
ylabel('Count');
xlim([0,100]);
subplot(1,3,2);
hist(r.^2*pi,0:0.3:10);
xlabel('PVS radius (Unit: Voxel)');
ylabel('Count');
xlim([0,7]);

subplot(1,3,3);hist(v,0:5:500);
xlim([0,200]);
set(gca,'FontSize',12);
xlabel('PVS volume (Unit: number of voxels)');
ylabel('Count');

%}




    
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
        
        
        
    

    
    



    
