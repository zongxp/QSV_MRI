function a=tof_bg_rm(dname,roi_name,mask_threshold,prefix,keep_constant,common_background)
% a=tof_bg_rm(dname,roi_name,mask_threshold,prefix,keep_constant,common_background)
% roi_name: name of the roi
% mask_threshold: used when roi_name is empty.  The threshold is calculated
% as max(dname(:))*mask_threshold
%
if ~exist('keep_constant','var')
    keep_constant=true;
end

if ~exist('common_background','var')
    common_background=false;
end

if isa(dname,'char')
    a=ri(dname);
else
    a=dname;
end

a=double(a);
if exist('roi_name','var')  &&~isempty(roi_name)
    if isa(roi_name,'char')
        m=ri(roi_name);
    else
        m=roi_name;
    end
else
    
    m=a>max(a(:))*mask_threshold;
    figure;imshow4(double(m),[],[1,size(m,3)]);
end

x=[];
y=[];
for sl=1:size(a,3)
    if sum(vec(m(:,:,sl)))==0
        continue;
    end
for k=1:size(a,4)
    
    
    if k==1
        pos=ind2subb(size(m(:,:,sl)),find(m(:,:,sl)>0));
        posx=pos(:,1);
        posy=pos(:,2);
        
        
        meanx=mean(posx);
        meany=mean(posy);
        
        rx=diff(min_max(posx));
        ry=diff(min_max(posy));
        
        
        posx=(posx-meanx)/rx;
        posy=(posy-meany)/ry;
        
        %   x=[pos,ones(size(pos,1),1),pos.^2,pos(:,1).*pos(:,2)];
        
        x=[posx,posy,ones(size(pos,1),1),posx.^2,posy.^2,posx.*posy];
        
    end
    if common_background
        if k==1
            a_tmp=mean(a(:,:,sl,:),4);
            y=double(a_tmp(m(:,:,sl)>0));
            b=x\y;
        end
    else
        a_tmp=a(:,:,sl,k);
        y=double(a_tmp(m(:,:,sl)>0));
        b=x\y;
    end

       
        %%{
        %figure;plot(y,x*b,'o');
        %hold on;plot([min(y),max(y)],[min(y),max(y)],'r-');
        %ylim([min(x*b),max(x*b)]);
        %}
        for i=1:size(a,1)
            for j=1:size(a,2)
               % a(i,j,sl,k)=a(i,j,sl,k)-[i,j,1,i^2,j^2,i.*j]./norm*b;
                 i2=(i-meanx)/rx;
                 j2=(j-meany)/ry;
                a(i,j,sl,k)=a(i,j,sl,k)-[i2,j2,1,i2^2,j2^2,i2*j2]*b;
            end
        end
        
        if keep_constant
        a(:,:,sl,k)=a(:,:,sl,k)+mean(y);
        end
        
    end
    
end

%%

if exist('prefix','var')  && ~isempty(prefix)
    save([prefix,'_detrend.mat'],'a');
end

