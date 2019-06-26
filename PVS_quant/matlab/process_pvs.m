function process_pvs(fname,T2name,skip_q)
% 6/10/2019: in PVSLength, delete ind (n*1); use (indc n*3) instead to
% speed up. delete path variable, save nvox_path; since
% path{i}=1:nvox_path(i)

if ~exist('skip_q','var')
    skip_q=false;
end

[dname,fname2]=fileparts(fname);

prefix=strtok(fname2,'.');

if ~exist(prefix,'dir')
    mkdir(prefix);
end
cd(prefix);
tpvs=tic;
prefix2=[filename(prefix),'_PVS'];

PVSLength(fullfile('..',filename(fname)),prefix2);  %now called ThinningPathFind(mfile)
toc(tpvs);

pathCurved_nofix(prefix2);  %fixedpath
toc(tpvs);

fixGap(prefix2);   %probv % 1. why no i_ind_path in the file? 2. why this function can run when using bsub?
PVS_orient(prefix2);
toc(tpvs);

if skip_q
    cd ..;
    return;
end
PVS_quantify(prefix2);
toc(tpvs);

%fprintf('pvs_vol_frac_ind(%s_stat.mat,%s)\n',prefix2,T2name);
pvs_vol_frac_ind(prefix2,T2name);
toc(tpvs);

fprintf('run mat2niigz(%s)\n',[prefix2,'.mat']);
mat2niigz([prefix2,'.mat'],'c',prefix2,true,fullfile('..',filename(fname)));
toc(tpvs);

cd('..');
%%

function PVSLength(fname,out_prefix)
% m: pvs mask
% prefix: file name of the output file

m=ri_d1(fname);
% skip_thinning=true;

% if ~skip_thinning
%     [c,c_rm,out_thin]=thining(m); %skip thinning, since it will produce angle pathes
% else
[c,sub_c]=clusterize2(m>0,2);

tic;
fprintf('find path for %d PVS\n',length(sub_c));
for i=1:length(sub_c)
    [nvox_path(i),sub{i}]=findPathMatlab(sub_c{i});
end

for i=1:length(sub)
    ind{i}=sub2indb(size(m),sub{i});
end

if 0 %for debug only
    pathroi=zeros(size(m));
    cpath=double(c>0);
    
    for i=1:length(nvox_path)
        
        pos_path=sub{i}(1:nvox_path(i),:);
        
        %a=sub2indb(size(m),pos_path);
        
        pathroi=set_val(pathroi,pos_path,i,size(m));
        cpath=set_val(cpath,pos_path,2,size(m));
        
    end
end
%prefix=strtok(prefix,'.');

c=int16(c);
orient=get_orient(fname);
save_pvs(orient,c,ind,sub,nvox_path,out_prefix);



function d=set_val(d,sub,val,sz)

a=sub2indb(sz,sub);
d(a)=val;



for i=1:size(pos)
    d(pos(i,1),pos(i,2),pos(i,3))=val;
end

function [nvox_path,pos,maxlen]=findPathMatlab(pos)

%% pos is the ind of voxels in the cluster.1*nvox
%% sz: matrix size.
%%{


%pos=ind2subb(sz,ind);

for i=1:3
[tmp,imax(i)]=max(pos(:,i));
[tmp,imin(i)]=min(pos(:,i));
end

iall=[imax,imin];
dist=0;
pair=[];
for i=1:6
    for j=i+1:6
        
        
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
  %  disp(ind);
end

sp(size(sp,1)+1:size(sp,2),:)=0;

sp=tril(sp + sp');
[maxlen,pathmax] = graphshortestpath(sp,pair(1),pair(2),'Directed',false);
nvox_path=length(pathmax);



xpath=setdiff(1:size(pos,1),pathmax);
pos=[pos(pathmax,:);pos(xpath,:)];



%{

function [npix_path,pos2,maxlen]=findPathMatlab(pos)

%% pos is the ind of voxels in the cluster.1*nvox
%% sz: matrix size.
%%{


%pos=ind2subb(sz,ind);


for i=1:3
    [tmp,imax(i)]=max(pos(:,i));
    [tmp,imin(i)]=min(pos(:,i));
end

%iall=unique([imax,imin]);
iall=([imax,imin]);
dist=0;
% pair=[];
% for i=1:length(iall)
%     for j=i+1:length(iall)
%         
%         disttmp=sos(pos(iall(i),:)-pos(iall(j),:));
%         if dist<disttmp
%             pair=[iall(i),iall(j)];
%             dist=disttmp;
%         end
%         
%     end
% end
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


xpath=setdiff(1:size(pos,1),pathmax);

pos2=[pos(pathmax,:);pos(xpath,:)];
npix_path=length(pathmax);

%}
function [p,w] = find_neighbors2(pos)

p=[];
w=[];
for i=1:size(pos,1)
    pos1=pos(i,:);
    for j=i+1:size(pos,1)
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
% remove path with >90 degree turns

a=load([prefix,'.mat']);
ind=a.ind;
sub=a.sub;

sz=size(a.c);
nfound=0;
nchecked=0;
nvox_path=a.nvox_path;
fixedpath=zeros(sz,'int16');
probclust=fixedpath;
niter=0;
while nchecked<length(nvox_path)
    
    niter=niter+1;
    if mod(nchecked+1,50)==0
        %    disp(nchecked);
    end
    
    i=nchecked+1;
    %%
    
    pospath=sub{i}(1:nvox_path(i),:);
    
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
            
            nvox_path(i)=[];
            ind(i)=[];
            sub(i)=[];
            
            
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

c=zeros(sz,'int16');
pathroi=zeros(sz);

for i=1:length(nvox_path)
    
    c(ind{i})=i;
    
    pathroi(ind{i}(1:nvox_path(i)))=i;
    
end


%save([prefix,'_Curv_debug.mat'],'debug.mat','pathroi','fixedpath','probclust');
orient=get_orient(prefix);
save_pvs(orient,c,ind,sub,nvox_path,prefix);



function save_pvs(orient,c,ind,sub,nvox_path,prefix)
orient.c=c;
orient.ind=ind;
orient.sub=sub;
orient.nvox_path=nvox_path;
save(prefix,'-struct','orient');

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

function fixGap(prefix)

% assume that path{i} is 1:length(path{i})
a=load([prefix,'.mat']);
nvox_path=a.nvox_path;

sz=size(a.c);
c=a.c;
ind=a.ind;
sub=a.sub;
%% if there is a gap between the voxel and its path, delete that voxel

probclust=0*c;
probv=zeros(size(c),'int16');

for i=1:length(nvox_path)
    
    while 1
        voxremoved=false;
        nvox=nvox_path(i);
        while  nvox<length(ind{i})
            
            dist=sos(repmat(sub{i}(nvox+1,:),[nvox_path(i),1])-sub{i}(1:nvox_path(i),:));
            if min(dist)>=2
                
                [tmp,indtmp]=min(dist);
                
                if hasgap(sub{i},nvox+1,indtmp)
                    if ~voxremoved
                        probclust(ind{i})=i;  %was 1
                    end
                    
                    probv(ind{i}(nvox+1))=i;
                    c(ind{i}(nvox+1))=0;
                    ind{i}(nvox+1)=[];
                    
                    sub{i}(nvox+1,:)=[];
                    voxremoved=true;
                    continue;
                end
            end
            nvox=nvox+1;
        end
        
        if ~voxremoved
            break;
        end
        %if voxremoved, check again.
    end
    
end

if (sum(probv(:))>0)
    
    fprintf('%d PVSs found that have branches\n',  length(unique(probv(:)))-1);
    fprintf('%d voxels removed to fix the problem\n',sum(probv(:)>0));
    
else
    disp('No braches found for any path.  Good!');
end
%
% c=zeros(sz,'int16');
%
% for i=1:length(nvox_path)
%     c(ind{i})=i;
% end

orient=get_orient(prefix);
save_pvs(orient,c,ind,sub,nvox_path,prefix);




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






function [omask,clusters]=clusterize2(mask,thr)
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

clusters=clusters(nvox_sort>=thr);

nclus = length(clusters);

for i=1:nclus
    
    for j=1:size(clusters{i},1)
        
        a=clusters{i}(j,:);
        omask(a(1),a(2),a(3)) =i;
        
        
        
    end
    
end
% sort clusters to be compatible with earlier version

    for j=1:length(clusters)
        ind=sub2indb(size(omask),clusters{j});
        
        [tmp,i_ind]=sort(ind);
        clusters{j}=clusters{j}(i_ind,:);
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













