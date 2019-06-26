function ThinningPathFind(mfile,do_clusterize)
% originally named PVSLength
% 
% 
% mfile: matfile name
% contains d: mask; voxsize: if not found, assume [1,1,1].

if ~exist('do_clusterize','var')
    do_clusterize= true;
end
    
try
    voxsize=ri(mfile,[],[],'voxsize');
catch
    warning('voxsize not found; set to [1,1,1]');
    voxsize=[1,1,1];
end
try
m=ri(mfile,[],[],'c');
catch
m=ri(mfile,[],[],'d');
    
end

if do_clusterize
    [c,c_rm,out_thin]=thining(m);
else
    [c,c_rm,out_thin]=thining_noclust(m);
end

% if max(tmp(:))>0
%     fprintf('%d clusters removed by the Skeleton3D algorithm\n',max(tmp(:)));
%     fprintf('Results saved in %s.raw',[fname,'RemovedClustersSkeleton3D']);
%     write_mhd(c_rm>0,[fname,'RemovedClustersSkeleton3D']);
% end

%c=bwlabeln(m>0);
npath=0;
i_ind_path={};
ind={};

for i=1:max(c(:))
    if mod(i,10)==0
        disp(i);
    end
    
    indi=find(out_thin==i);
    if ~isempty(indi)
        npath=npath+1;
     [path_ind,ind_tmp]=findPathMatlab(indi,size(m));   
     
     tmp=find(c==i);
     
     c(tmp)=npath;
     
    tmp2=ind_tmp(path_ind);
    
    xind=setdiff(tmp,tmp2);
    
    ind{npath}=[tmp2;xind];
    i_ind_path{npath}=1:length(path_ind);
    
    end
    
end

pathroi=zeros(size(m));
cpath=double(c>0);
for i=1:length(i_ind_path)
    pathroi(ind{i}(i_ind_path{i}))=i;
    cpath(ind{i}(i_ind_path{i}))=2;
end

mfile=strtok(mfile,'.');
save(['ThPth_',mfile],'pathroi','c','ind','i_ind_path','cpath','c_rm','voxsize');
% c_rm: removed clusters by thinnig algorithm
% ind: 1*npath cell array containing the index for all voxels in a cluster
% i_ind_path: 1*npath cell array, containing the indices to pick ind elements for the path voxels
% cpath: mask with path voxels set to 2 and other voxels with c>0 set to 1.
% c: clusters after thinning algorithm
% pathroi: mask for voxels in the path

 
function [out,out_rm,out_thin]=thining(c) 

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
       

function [out,out_rm,out_thin]=thining_noclust(c) 

    out=c;
    
out_rm=0*out;
out_thin=0*out;

for i=1:max(out(:))
   
    if mod(i,50)==0
        disp(i);
    end
   
    p=Skeleton3D(out==i);
    
    m=(out==i)&p>0;
    mi=(out==i);
    if sum(m(:))==0
      out(mi)=0;
      out_rm(mi)=1;
      
    else
      
        out_thin(m)=i;
    end
    
end

function write_mhd(d,fname)

roi=permute(d,[2,1,3]);

f=fopen([fname,'.raw'],'wb');
%f=fopen(fname,'rb');

d=fwrite(f,roi,'integer*2');

fclose(f);
          

function [pathmax,ind,maxlen]=findPathMatlab(ind, sz)

%% ind is the ind of voxels in the cluster.1*nvox
%% sz: matrix size.
%pathmax: the indices for the path voxels
% ind: the indices for the 



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

sp=sparse(s(:,1)',s(:,2)',w);  % max of s(:,1) <= max of s(:,2)
sp(size(sp,1)+1:size(sp,2),:)=0;

sp=tril(sp + sp');  % what is purpose for this?

[maxlen,pathmax] = graphshortestpath(sp,pair(1),pair(2),'Directed',false);

ind=ind(:);


xpath=setdiff(1:length(ind),pathmax);
ind=[ind(pathmax);ind(xpath)];
pathmax=1:length(pathmax);


    function [p,w] = find_neighbors2(pos)
% find all neighbor pairs and set their weights to w
       
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
    
    
    
        
        
    