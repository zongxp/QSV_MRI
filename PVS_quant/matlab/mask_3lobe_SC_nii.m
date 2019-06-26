function mask_3lobe_SC_nii(fnii)

%freesurfer template
% the data set should be lr (1 dim); ant to pos (2nd); inf to sup (3rd);
% 
% if isa(mat_fname,'char') 
%     a=ri(mat_fname);
% else
%     a=mat_fname;
% end

nii=load_untouch_nii(fnii);

a=nii.img;
a=permute(a,[1,3,2]);  % make the second dimension ant to pos

%a=permute(a,[2,1,3]);

m_sc=mask_subcortical(a);
m_frn=(a>=1&a<=20) | (a>=23 & a<=26);
m_par=(a>=57 & a<=62);

m_oc=a>=43 & a <=54;
m_tempr(:,:,:,1)=double(a==79 | a==81| a==83 | a==85 | a==87 | a==89  |a==55);
m_tempr(:,:,:,2)=double(a==80 | a==82| a==84 | a==86 | a==88 | a==90  |a==56);

m_insu = a==29|a==30;
m_cing=a>=31& a<=36;
m_hip=a>=37 & a<=40;
m_postL=a==57;
m_postR=a==58;
m_post=m_postL|m_postR;

m_precL=a==1; 
m_precR= a==2;
m_prec=m_precL|m_precR;

%% for temporal lobe fill the holes

for i=1:size(m_tempr,1)  
    for j=size(m_tempr,2)-1:-1:1 %  pos to ant
       disp([i,j]);
        for k=size(m_tempr,3)-1:-1:1  % sup to inf
            
            for l=1:2
              if m_tempr(i,j,k,l)==0 && m_tempr(i,j,k+1,l)>0 && any(m_tempr(i,j+1:end,k,l)~=0)
                m_tempr(i,j,k,l)=1;
              end
            end
        end
    end
    
end

m_tempr=double(m_tempr(:,:,:,1))+double(m_tempr(:,:,:,2))*2;

brainmask=a>=1&a<=90;

for i=1:size(brainmask,2)
    for j=1:size(brainmask,3)
        
        ind=find(brainmask(:,i,j)>0);
        if isempty(ind)
            continue;
        end
        brainmask(ind(1):ind(end),i,j)=1;
    end
end

        
%writeraw(uint8(brainmask),'mask_brain.raw');





%m=m_sc*1+m_frn*2+m_par*3+m_tempr*4+m_oc*5+m_insu*6+m_cing*7+m_hip*8;
%writeraw(uint8(m),'mask_roi.raw');

m2=m_prec+m_post*2;
%writeraw(uint8(m2),'mask_pre_post_centr.raw');

res=m2*0;


for i=size(a,3):-1:1
   
   indyL=ind2subb(size(a(:,:,i)),find(m_precL(:,:,i)>0));
   indyR=ind2subb(size(a(:,:,i)),find(m_precR(:,:,i)>0));

   if isempty(indyL) || isempty(indyR)
       if ~exist('last_thr','var')
           continue;
       end
         thr=last_thr;
   else
       
       
   y_maxL=max(indyL(:,2));
   y_maxR=max(indyR(:,2));
   
   thr=round((y_maxL+y_maxR)*0.5);
   last_thr=thr;   
 %  disp(last_thr);
   end
   
   tmp=a(:,:,1)*0;
   tmp(:,1:thr)=1;  % frontal
   tmp(:,thr+1:end)=2; % pariental occ
   
   tmp(m_prec(:,:,i)>0)=1;
   tmp(m_post(:,:,i)>0)=2;
   tmp(m_tempr(:,:,i)>0)=3;
   
   res(:,:,i)=tmp;
       
end
res(m_sc>0)=4;
res=res.*brainmask;

%writeraw(uint8(res),'mask_lobes.raw');
nii.img=res;
prefix=strtok(fnii,'.');
save_untouch_nii(nii,[prefix,'_3lobe_SC']);
%save mask_lobes res
%save mask_temporal_lobe m_tempr





function m3=mask_subcortical(d)

m=(d>=9&d<=13) | (d>=48&d<=52);
m2=bwmorph3d(m,'dilate',15);

m3=bwmorph3d(m2,'shrink',10);
m3=m3|m;

