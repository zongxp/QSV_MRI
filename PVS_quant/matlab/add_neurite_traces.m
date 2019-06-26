function [V,fname]=add_neurite_traces(data_fname,dir_trace,xform,sig_scale,pathfile,xform_out)
% xform: the transformation needed to match data_fname to imageJ image.
% sig_scale: if <0; then assume dark lines.  default 2.
% if sig_scale is empty, then no extra voxels will be added to the mask,
% only the voxels that intersect with the path will be included in the ROI.
% the V is correct only if their is no U-turns in the path defined by
% neurite tracer.

if ~exist('sig_scale','var')
    sig_scale=2;
end

if ~exist('xform_out','var')
    xform_out={'xy','flipy'};
end
%always use dicom images.  dicom and .mat files matches in matrix
    % size, while .hdr files may not.
    [d,dm]=ri_mat_dcm_ana(data_fname,[]);
    
  %  dm=ri(data_fname,[],[],'voxsize');
   
%d=dataxform(d,xform,0);
dm= dimxform(dm,xform,0);

sz=size(d);
roi=zeros(sz,'single');
roi_path=zeros(sz,'single');

m_background=zeros(sz,'single');

V = zeros([sz,3],'single');

tort = zeros(sz,'single');

dir_str=dir(fullfile(dir_trace,'*.swc'));

path=cell(1,length(dir_str));
ind=cell(1,length(dir_str));

orient = cell(1,length(dir_str)); % n*3 vectors; the two end points on each side are 0 since it is not well defined.
tort_seg = cell(1,length(dir_str));
curv_seg = cell(1,length(dir_str));
tort_path=zeros(1,length(dir_str));
fprintf('adding %d paths\n',length(dir_str));

for j=1:length(dir_str)
    fname=fullfile(dir_trace,dir_str(j).name);
    fid=fopen(fname,'r');
    
    
    a=textscan(fid,'%f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);
    pos=zeros(length(a{1}),3);
    for i=1:length(a{1})
        pos(i,:)=[a{3}(i),a{4}(i),a{5}(i)]./dm+1;
    end
    
    pos=posxform(pos,xform,sz,1);
    
    
    pos=round(connect_pos(pos));
    
    orient{j} = zeros(size(pos));
    
    tort_seg{j} = zeros(1,size(pos,1));
    curv_seg{j} = zeros(1,size(pos,1));
    
    
    for i=1:size(pos,1)
        
        if i<3 || i>size(pos,1)-2
            continue;
        end
        
        pca_res=pca(pos(i-2:i+2,:));
        orient{j}(i,:) = pca_res(:,1);
        
        tort_seg{j}(i) = tortuosity(pos(i-2:i+2,:),dm);
        
        if i<5 || i>size(pos,1)-4
            continue;
        end
        
        curv_seg{j}(i) = curvature(pos(i-4:i+4,:),dm);
        
    
    end
    tort_path(j) = tortuosity(pos,dm);
    %% find the direction vector for the path
    
    
    ind{j}=sub2indb(sz,pos);
    path{j}=1:size(pos,1);  % the initial voxels are path voxels
    
    
    pos=pos.*repmat(dm,[size(pos,1),1]);
    
    %%{
    if ~isempty(sig_scale)
        me=neighbors_mask(pos,1.2,sz,dm);
        men=neighbors_mask(pos,0.8,sz,dm);
        
        
        
        me=me&~men;
        
        mn=mean(d(me));
        sd=std(d(me));
        
        if sig_scale>0
            roi(d>(mn+sig_scale*sd) & men)=j;
        else
            roi(d<(mn+sig_scale*sd) & men)=j;
        end
    end
    %}
    for i=1:size(pos,1)
        id=round(pos(i,:)./dm);
        roi(id(1),id(2),id(3))=j;
        roi_path(id(1),id(2),id(3))=j;
        
        m_background(id(1),id(2),id(3))=2;
        
        V(id(1),id(2),id(3),:)=orient{j}(i,:);
        tort(id(1),id(2),id(3))=tort_seg{j}(i);
         
    end
    
    if ~isempty(sig_scale)
        roi_extra=(roi==j)&(roi_path~=j);
        ind_extra=find(roi_extra(:)>0);
        ind{j}(end+1:end+length(ind_extra))=ind_extra;
        
        %m=m|(d>(mn+2*sd) & men);
        m_background(me>0)=1+m_background(me>0);
    end
end


% if ~isempty(strfind(data_fname,'.hdr'))
%     roi=flip(permute(roi,[2,1,3]),2);
%     roi_path=flip(permute(roi_path,[2,1,3]),2);
% end

roi=dataxform(roi,xform_out,0);
roi_path=dataxform(roi_path,xform_out,0);


c=single(roi);
roi_path=single(roi_path);
%save(sprintf('PVS_%s_background',prefix_trace),'m_background');
voxsize=dm;
dm2=reshape(dm,[1,1,1,3]);
V2 = V.*repmat(dm2,[sz,1]);

 phi=atan2(V2(:,:,:,2),V2(:,:,:,1))*180/pi;
 th=acos(V2(:,:,:,3)./sqrt(sos(V2,4)))*180/pi;
    
phi=single(phi);
th=single(th);
phi(isnan(phi))=0;
th(isnan(th))=0;

V = single(V);

if isempty(sig_scale)
    fname=pathfile;
save(fname,'ind','c','path','roi_path','voxsize','orient','phi','th','V','tort','tort_seg','tort_path','curv_seg');
else
    fname=sprintf('%s_sig%s.mat',pathfile,num2str(sig_scale));
    save(fname,'ind','c','path','roi_path','voxsize','orient','phi','th','V');
end
%m=flipdim(m,2);

%roi(m>0)=m(m>0);
function m=neighbors_mask(pos,dist,sz,dim)

m=zeros(sz);
for i=1:size(pos,1)
    
    for j1=-ceil(dist/dim(1)):ceil(dist/dim(1))
        for j2=-ceil(dist/dim(2)):ceil(dist/dim(2))
            for j3=-ceil(dist/dim(3)):ceil(dist/dim(3))
                
                
                ind=round(pos(i,:)./dim)+[j1,j2,j3];
                
                if any(ind<1) || any(ind>sz)
                    continue;
                end
                
                if m(ind(1),ind(2),ind(3))==1
                    continue;
                end
                
                if sos([j1,j2,j3].*dim)<=dist
                    m(ind(1),ind(2),ind(3))=1;
                end
            end
        end
    end
end

%{
old implementation; has some issues
function m=neighbors_mask(pos,dist,sz,dim)
%pos are one-based integers; n*3


pos2=pos(1,:);  %pos2 removes duplicate entries in pos.

for i=2:size(pos,1)
    if ~any(pos(i,:)~=pos2(end,:))
        continue;
    end
    pos2(end+1,:)=pos(i,:);
end

pos=pos2;

th=(0:19)*2*pi/20;

dpos=diff(pos,1,1);

m=zeros(sz);
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
    
 orig=pos(i:i+1,:);
 for k=1:length(th)
    
     for l=1:size(orig,1)
      coord=round((orig(l,:)+dist*(n2*sin(th(k))+n3*cos(th(k))))./dim);
      
      coord(coord<1)=1;
      coord(coord>sz)=sz(coord>sz);
      
      m(coord(1),coord(2),coord(3))=1;
      
     end
 end
 
 
end
    
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



