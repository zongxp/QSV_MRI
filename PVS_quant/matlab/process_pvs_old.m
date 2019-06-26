function process_pvs_old(fname,T2name)
  
    [dname,fname2]=fileparts(fname);
    
    prefix=strtok(fname2,'.');
    tpvs=tic;
    if ~exist(fullfile(dname,prefix),'dir')
       mkdir(fullfile(dname,prefix));
    end
       cur_d= cd(fullfile(dname,prefix));
            
        prefix2=[filename(prefix),'_PL'];
        if ~exist([prefix2,'.mat'],'file')
         disp(fname);
         PVSLength(fullfile('..',filename(fname)),prefix2);  %now called ThinningPathFind(mfile)  
         toc(tpvs);
        end
        
        if ~exist([prefix2,'_Curv.mat'],'file')
         disp(prefix2);
         pathCurved_nofix(prefix2);  %fixedpath
           toc(tpvs);
        end
        
        if ~exist([prefix2,'_Curv_Dv.mat'],'file')
         disp([prefix2,'_Curv.mat']);
         fixGap([prefix2,'_Curv']);   %probv % 1. why no i_ind_path in the file? 2. why this function can run when using bsub?
           toc(tpvs);
        end
        
        if ~exist([prefix2,'_Curv_c.nii.gz'],'file')
          mat2niigz([prefix2,'_Curv.mat'],'c',[prefix2,'_Curv_c'],true,fullfile('..',filename(fname)));
            toc(tpvs);
        end
        
        if ~exist([prefix2,'_Curv_Dv_stat.mat'],'file')
            disp([prefix2,'_Curv_Dv.mat']);
            PVS_quantify([prefix2,'_Curv_Dv.mat']);      
              toc(tpvs);
        end
        
       if ~exist([prefix2,'_Curv_Dv_stat_orient.mat'],'file')
           disp([prefix2,'_Curv_Dv_stat.mat']);
         PVS_orient([prefix2,'_Curv_Dv_stat.mat']);
           toc(tpvs);
       end
       
       if ~exist([prefix2,'_Curv_Dv_stat_vfind.mat'],'file') && exist('T2name','var')
           
           fprintf('pvs_vol_frac_ind(%s_Curv_Dv_stat.mat,%s)\n',prefix2,T2name);
           pvs_vol_frac_ind([prefix2,'_Curv_Dv_stat.mat'],T2name);
             toc(tpvs);
           
       end
       
    cd(cur_d);
    


%%

function PVSLength(fname,out_prefix)
% m: pvs mask
% prefix: file name of the output file

m=ri_d1(fname);
skip_thinning=true;

if ~skip_thinning
    [c,c_rm,out_thin]=thining(m); %skip thinning, since it will produce angle pathes
else
    c=clusterize2(m>0,2);
    out_thin = c;
end
%{
if max(tmp(:))>0
    fprintf('%d clusters removed by the Skeleton3D algorithm\n',max(tmp(:)));
    fprintf('Results saved in %s.raw',[prefix,'_RemovedClustersSkeleton3D']);
    write_mhd(c_rm>0,[prefix,'_RemovedClustersSkeleton3D']);
    drawnow;
end
%}
%c=bwlabeln(m>0);
npath=0;
path={};
ind={};
for i=1:max(c(:))
    if mod(i,10)==0
        fprintf('%d/%d\n',i,max(c(:)));
    end
    
    indi=find(out_thin==i);
    if ~isempty(indi)  && length(indi)>1
        npath=npath+1;
     [path_ind,ind_tmp]=findPathMatlab(indi,size(m));   
     
     tmp=find(c==i);
     
     c(tmp)=npath;
     
    tmp2=ind_tmp(path_ind);
    
    xind=setdiff(tmp,tmp2);
    
    ind{npath}=[tmp2;xind];
    path{npath}=1:length(path_ind);
    
  
    end
    
end

pathroi=zeros(size(m));
cpath=double(c>0);
for i=1:length(path)
    pathroi(ind{i}(path{i}))=i;
    cpath(ind{i}(path{i}))=2;
end

%prefix=strtok(prefix,'.');

save([out_prefix,'.mat'],'pathroi','c','ind','path','cpath');

copy_orient(fname,[out_prefix,'.mat']);
 
function [out,out_rm,out_thin]=thining(c) 

    % out: clusterizes, excluding those removed by the thinning algorithm
    % out_rm: clusters removed by the thinning algorithm
    % out_thin: thinned clusters
    
    out=clusterize2(c>0,2);
    p=Skeleton3D(c>0);

out_rm=0*out;
out_thin=0*out;

for i=1:max(out(:))
   
    if mod(i,50)==0
        disp(i);
    end
    
    m=(out==i)&p>0;
    mi=(out==i);
    if sum(m(:))==0
      out(mi)=0;
      out_rm(mi)=1;
      
    else
      
        out_thin(m)=i;
    end
    
end
       



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

try
sp=sparse(s(:,1)',s(:,2)',w);

catch
    disp('s');
    disp(s);
    disp('ind');
    disp(ind);
end

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
    
    
function res = sos(x ,dim, pnorm)
% res = sos(x [,dim, pnorm])
%
% function computes the square root of sum of squares along dimension dim.
% If dim is not specified, it computes it along the last dimension.
%
% (c) Michael Lustig 2009

if nargin < 2
    dim = size(size(x),2);
end

if nargin < 3
    pnorm = 2;
end


res = (sum(abs(x.^pnorm),dim)).^(1/pnorm);
    
 function pathCurved_nofix(prefix)
% file should contain a ind file (1*npath)
% path file (1*npath)
% c: non zero with values equal to the npath index.

a=load([prefix,'.mat']);
ind=a.ind;
sz=size(a.c);
nfound=0;
nchecked=0;
path=a.path;
fixedpath=zeros(sz);
probclust=fixedpath;
niter=0;
while nchecked<length(path)
   
    niter=niter+1;
    if mod(nchecked+1,50)==0
    %    disp(nchecked);
    end
    
    i=nchecked+1;
   %% 
    xpath=setdiff(1:length(ind{i}),path{i});
     
       ind{i}=ind{i}(:);
    ind{i}=[ind{i}(path{i});ind{i}(xpath)];
 
    path{i}=1:length(path{i});
    
   pospath=ind2subb(sz,ind{i}(path{i}));
    
   breakpath=false;
   for j=1:size(pospath,1)-2
       
       d1=pospath(j,:)-pospath(j+1,:);
       d2=pospath(j+1,:)-pospath(j+2,:);  %was j+2 and j+3 before; fixed 12/19/2018

   
       cs=sum(d1.*d2)/sos(d1)/sos(d2);
       
       if cs<0
           probclust(ind{i})=1;
           
           probclust(pospath(j,1),pospath(j,2),pospath(j,3))=2;
           probclust(pospath(j+1,1),pospath(j+1,2),pospath(j+1,3))=2;
           
           
           probclust(pospath(j+2,1),pospath(j+2,2),pospath(j+2,3))=2;
       
           path(i)=[];
           ind(i)=[];
           
    
           
           nfound=nfound+1;
           breakpath=true;
           break;
           
       end
       
   end
 %%  
   if breakpath
       continue;
   end
   
   
   
   nchecked=nchecked+1;
   
end

if nfound>0
fprintf('%d paths found with >90 degree angle\n',nfound);

else
    disp('No paths with > 90 degree angle found.  Good!');
end

c=zeros(sz);
pathroi=zeros(sz);

for i=1:length(path)

    c(ind{i})=i;
    
    pathroi(ind{i}(path{i}))=i;

end


save([prefix,'_Curv.mat'],'c','ind','pathroi','path','fixedpath','probclust');

  copy_orient(prefix,[prefix,'_Curv.mat']);
    

function sub=ind2subb(mtrx,ind_arr)

sub=zeros(length(ind_arr),length(mtrx));
for j=1:length(ind_arr)
     ind=ind_arr(j);
for i=1:length(mtrx)
   
    if i==length(mtrx)
        sub(j,1)=ind;
    else
    sub(j,end-i+1)=floor((ind-1)/prod(mtrx(1:end-i)))+1;
    end
    ind=ind-(sub(j,end-i+1)-1)*prod(mtrx(1:end-i));
end

end

    function res=minDist(ind0,pathind,sz)
        
          pos1=ind2subb(sz,ind0);
          pos2=ind2subb(sz,pathind);
            
          dist=sos(repmat(pos1,[size(pos2,1),1])-pos2);
            
          res=min(dist);
                
                
        


function fixGap(prefix)

% assume that path{i} is 1:length(path{i})
a=load([prefix,'.mat']);
path=a.path;

sz=size(a.c);
c=a.c;
ind=a.ind;
%% if there is a gap between the voxel and its path, delete that voxel

probclust=0*c;

nspath=0;
rpath=[];
for i=1:length(path)
  if length(path{i})==1
      rpath(end+1)=i;
      fprintf('Path %d has only 1 voxel\n',i);
      probclust(ind{i})=1;
      nspath=nspath+1;
  end
    
end

if nspath>0
path(rpath)=[];
ind(rpath)=[];
fprintf('%d paths has only 1 voxel\n',nspath);
end
%%
probv=zeros(size(c));
for i=1:length(path)
    
    if mod(i,50)==0
     %   disp(i);
    end
    pos=ind2subb(sz,ind{i});
    
    pospath=ind2subb(sz,ind{i}(path{i}));
        
    while 1
        voxremoved=false;
      nvox=length(path{i});
      while  nvox<length(ind{i})
            
            dist=sos(repmat(pos(nvox+1,:),[size(pospath,1),1])-pospath);
            
            
      
            if min(dist)>=2
                
                [tmp,indtmp]=min(dist);
                
                if hasgap(pos,nvox+1,path{i}(indtmp))
                  
                  if ~voxremoved  
                    probclust(ind{i})=i;  %was 1
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

  fprintf('%d PVSs found that have branches\n',  length(unique(probv(:)))-1);
  fprintf('%d voxels will be removed to fix the problem\n',sum(probv(:)>0));
  %fprintf('%d voxels will be deleted from %d paths; paths:\n',sum(probv(:)>0),length(unique(probv(:)))-1);

else
    disp('No braches found for any path.  Good!');
end

%disp(unique(probv(probv>0))');



c=zeros(sz);
pathroi=zeros(sz);

for i=1:length(path)

    c(ind{i})=i;
    
    pathroi(ind{i}(path{i}))=i;

end

save([prefix,'_Dv.mat'],'ind','c','path','probv','pathroi','probclust');

 copy_orient(prefix,[prefix,'_Dv.mat']);



    
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
    
    function ind=pairIndex(i,j)
        
        if (i>j)
            tmp=j;
            j=i;
            i=tmp;
        end
        ind=i*10000+j;
        
        
  function omask=clusterize2(mask,thr)
%[mask,clusters]=clusterize2(mask,thr)
% find connected clusters in the mask and remove voxels with cluster size
% less than thr.
% only version; NO clusterize or clusterize1.


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
        
          nb = find_neighbors_clust(cluster_temp(ivc,:),mask);
        
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
    
        
        
function nb = find_neighbors_clust(pos,mask)

nb=[];
   sz = size(mask);
   if numel(sz) <3
       sz(3)=1;
   end
  
  for i=-1:1
      for j=-1:1
          for k=-1:1
            if(abs(i)==1 &&abs(j)==1 && abs(k)==1)
              %  continue;
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
      

    
   
     
   
      
    
   
    
        
        
    
