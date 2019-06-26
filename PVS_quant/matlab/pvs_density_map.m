function pvs_density_map(fname,dia,fwm)
% max_dist: diameter of the sphere for calculating the mean diameter

c=ri(fname,'','','c');
ind=ri(fname,'','','ind');
path=ri(fname,'','','path');
dpath=ri(fname,'','','dpath');

mwm=ri(fwm);

prefix=strtok(fname,'.');
%cpath=get_cpath(c,path,ind);  % path mask; value = path index
dmap=get_diam_map(c,path,ind,dpath);                              % diameter map;

%% 

%msmall=mask_sphere(dia+1,dia/2,[dia/2+1,dia/2+1,dia/2+1],true);

msmall=mask_sphere_weighted(dia+1,dia/2,[dia/2+1,dia/2+1,dia/2+1],true);


%%

res=c*0;

n=c*0;

for i=1:size(c,1)

    tic;
    for j=1:size(c,2)
        for k=1:size(c,3)
           
        if dmap(i,j,k)==0
            continue;
        end
         
        if abs(i-165)>=dia/2 || abs(j-212)>=dia/2 || abs(k-116)>=dia/2
          % continue; 
        end
        [res,n]=setv_img(res,[i,j,k],msmall,dmap(i,j,k),n);    
        
           %m=get_mask(size(c),[i,j,k],dia,msmall);
%          [c_crp,mcrop]=crop_img(c,[i,j,k],dia,msmall);
%         
%          mroi=mcrop&c_crp>0;
%          if sum(mroi(:))==0
%              continue;
%          end
         
          % dmap_crp=crop_img(dmap,[i,j,k],dia);
          % res(i,j,k)=mean_roi(dmap_crp,mroi);  
         
        end
    end
    fprintf('%d/%d: time remaining: %5.1f s\n',i,size(c,1),(size(c,1)-i)*toc);
end
n2=n.*double(mwm>0);
res2=0*res;
res2(n2>0)=res(n2>0)./n2(n2>0);

%res2=res2.*double(mwm>0);
save(sprintf('%s_diamMap_d%s.mat',prefix,num2str(dia)),'res2','n2');

function   [c,n]=setv_img(c,cen,msmall,val,n)

dia=size(msmall)-1;
sz=size(c);
i1=cen(1)-dia(1)/2:cen(1)+dia(1)/2;
i2=cen(2)-dia(2)/2:cen(2)+dia(2)/2;
i3=cen(3)-dia(3)/2:cen(3)+dia(3)/2;

j1=i1>=1&i1<=sz(1);
j2=i2>=1&i2<=sz(2);
j3=i3>=1&i3<=sz(3);

i1=i1(j1);
i2=i2(j2);
i3=i3(j3);


c(i1,i2,i3)=c(i1,i2,i3)+val*msmall(j1,j2,j3);

n(i1,i2,i3)=n(i1,i2,i3)+msmall(j1,j2,j3);

function m=mask_sphere_weighted(dim,rad,center,include_equal)
% m=mask_circle(dim,rad,center,include_equal)
if ~exist('include_equal','var')
    include_equal=true;
end

if length(dim)==1
    dim=dim*ones(1,3);
end

[x,y,z]=meshgrid(1:dim(2),1:dim(1),1:dim(3));

r2=(x-center(2)).^2+(y-center(1)).^2+(z-center(3)).^2;
if include_equal
    m=(r2<=rad^2).*(1+cos(pi*sqrt(r2)/rad))/2;
else
    m=(r2<rad^2).*(1+cos(pi*sqrt(r2)/rad))/2;
end



function [res,mcrop]=crop_img(m,cen,dia,msmall)

sz=size(m);
i1=cen(1)-dia/2:cen(1)+dia/2;
i2=cen(2)-dia/2:cen(2)+dia/2;
i3=cen(3)-dia/2:cen(3)+dia/2;

j1=i1>=1&i1<=sz(1);
j2=i2>=1&i2<=sz(2);
j3=i3>=1&i3<=sz(3);

i1=i1(j1);
i2=i2(j2);
i3=i3(j3);


res=m(i1,i2,i3);
if exist('msmall','var')
mcrop=msmall(j1,j2,j3);
end

function m=get_mask(sz,cen,dia,msmall)


i1=cen(1)-dia/2:cen(1)+dia/2;
i2=cen(2)-dia/2:cen(2)+dia/2;
i3=cen(3)-dia/2:cen(3)+dia/2;

j1=i1>=1&i1<=sz(1);
j2=i2>=1&i2<=sz(2);
j3=i3>=1&i3<=sz(3);

i1=i1(j1);
i2=i2(j2);
i3=i3(j3);

m=zeros(sz);
m(i1,i2,i3)=msmall(j1,j2,j3);




function res=get_diam_map(c,path,ind,dpath)

res=c*0;
for i=1:length(ind)
  res(ind{i}(path{i}))=dpath{i};    
end



function res=get_cpath(c,path,ind)

res=c*0;

for i=1:length(ind)
   res(ind{i}(path{i}))=i;    
end

