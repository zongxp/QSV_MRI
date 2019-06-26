function PVS_quantify(fname)
%fname = 'PVSLength_testzb_Curv_Dv_Comb_Dv.mat';
%fname='PVSLength_HVB238b_Curv_Dv_Comb_Dv.mat';
%fname='PVSLength_HVB239c_Curv_Dv_Comb_Dv.mat';
% vsize: 1*3 array; voxel size

fprintf('Run PVS_quantify(%s) ...\n',fname);
if strcmp(fname(end-3:end),'.mat')
    fname=fname(1:end-4);
end

a=load(fname);
ind=a.ind;
voxsize=a.voxsize;

nvox_path=a.nvox_path;
npath=length(nvox_path);

ind_pathVoxelLabel = cell(1,npath); % which path voxel does each element in ind associate with (min distance)

c=a.c;
sub=a.sub;

l=zeros(1,npath);
lnorm=cell(1,npath);
v=zeros(1,npath);
r=zeros(1,npath);

v2=[];  % v for segments

dpath=cell(1,npath);
vpath=cell(1,npath);
%lpath=cell(1,npath);
dpath_all=[];
lthr=0.8; %lthreshold = 3 mm
diam=c*0;

for i=1:npath
 
    posc=sub{i};  
    posc=posc.*repmat(voxsize(:)',[size(posc,1),1]);
    
    orient=pca(posc);
   
 
    pos=posc(1:nvox_path(i),:);


    v(i)=length(ind{i})*prod(voxsize);
    lnorm{i}=zeros(1,nvox_path(i));
    
    for j=1:nvox_path(i)-1
        l(i)=l(i)+sqrt(sum((pos(j,:)-pos(j+1,:)).^2));
        lnorm{i}(j+1)=l(i);   
    end
    
    lnorm{i}=lnorm{i}/l(i);
    
        l(i)=l(i)+sos(pos(1,:)-pos(2,:))/2+sos(pos(end-1,:)-pos(end,:))/2;
        r(i)=sqrt(v(i)/l(i)/pi);
    
      pospath=pos;
        
      rpath_tmp=zeros(1,nvox_path(i)-1);
      rpath{i}=zeros(1,nvox_path(i));
      

      lseg=zeros(1,nvox_path(i));
      for j=2:nvox_path(i)-1
         
          lseg(j)=0.5*(sos(pospath(j,:)-pospath(j-1,:)))+0.5*(sos(pospath(j,:)-pospath(j+1,:)));
          
      end
      lseg(1)=sos(pospath(2,:)-pospath(1,:));
      lseg(end)=sos(pospath(end,:)-pospath(end-1,:));
      
      vseg=zeros(1,nvox_path(i));
      vsegclust=zeros(1,nvox_path(i));
      ind_pathVoxelLabel{i} = zeros(1,nvox_path(i));
      
      for j=1:length(ind{i})
            
          dist=sos(repmat(posc(j,:),[size(pos,1),1])-pos);
          
          [mdist,ind_min]=min(dist);
          vseg(mdist==dist)=vseg(mdist==dist)+prod(voxsize)/sum(mdist==dist);
          
          ind_pathVoxelLabel{i}(j)=ind_min;
         % if ind_min>2&&ind_min<size(pos,1)-1
         try
           ctheta_seg(ind{i}(j)) = theta_seg{i}(ind_min);
           cphi_seg(ind{i}(j)) = phi_seg{i}(ind_min);  
         catch
            disp(''); 
         end
         % end
      end
      
      dpath{i}=sqrt(vseg./lseg/pi)*2;
      
       for j=1:length(ind{i})
          dist=sos(repmat(posc(j,:),[size(pos,1),1])-pos);         
          [mdist,ind_min]=min(dist);
          vsegclust(j)=mean(vseg(mdist==dist));
          diam(ind{i}(j))=mean(dpath{i}(mdist==dist));
       end
      
      v2(end+1:end+length(vseg))=vseg;
       
      
      vpath{i}=vseg;
      vclust{i}=vsegclust;
    
      if l(i)>=lthr
       dpath_all=[dpath_all,dpath{i}]; % combine them so that we can plot below
      end
      
  
end

% only l, v , and dpath are used in pvs_stat

save([fname,'_stat'],'l','v','dpath','lthr','dpath_all','r',...
    'nvox_path','ind','vpath','v2','vclust','lnorm','ind_pathVoxelLabel','diam');


function ind=neighbors_mask(pos,dist,sz)
%pos are one-based integers; n*3

th=(0:19)*2*pi/20;

dpos=diff(pos,1,1);
m=zeros(sz);
ind=[];
for i=1:size(pos,1)-1
   n1=dpos(i,:);
 [tmp,t]=max(abs(n1));
 
 j=1:3;
 j(t)=[];
 n2=ones(1,3);
  
 n2(t)=-sum(n1(j).*[1,1])/n1(t);
 
 n2=n2/sqrt(sum(n2.^2)); 
 n1=n1/sqrt(sum(n1.^2));
 n3=cross(n1,n2);
    
 orig=zeros(5,3);
 ds=linspace(0.01,0.99,5);
 for k=1:5
 orig(k,:)=pos(i,:)+dpos(i,:)*ds(k);  % avoid passing the corners 
end
 
 for k=1:length(th)
    for j=0:ceil(dist)
     for l=1:size(orig,1)
      coord=round(orig(l,:)+j*(n2*sin(th(k))+n3*cos(th(k))));
      
      coord(coord<1)=1;
      coord(coord>sz)=sz(coord>sz);
      
  %    m(coord(1),coord(2),coord(3))=1;
      ind(end+1)=sub2indb(sz,coord);
     end
    end
 end
 
 
end
    
ind=unique(ind);

    
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


function check_dir(pos)

for i=1:size(pos,1)-2
   
 n=pos(i+1,:)-pos(i,:);
 n2=pos(i+2,:)-pos(i+1,:);
 
 c=sum(n.*n2)/sos(n)/sos(n2);
%disp([i,c]);
 if c<0
     disp([i,c]);
     disp('path error');
 end
 
 
end

function [Im]=sos(iX)

% room sum square of the last non-single dim.

if ndims(iX)>2

   Im=sqrt(sum((abs(iX).^2),ndims));

else
    
    if size(iX,2)==1
     Im=sqrt(sum((abs(iX).^2),1));
    else
    Im=sqrt(sum((abs(iX).^2),2));
    end
end

