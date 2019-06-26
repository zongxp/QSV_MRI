function pvs_colorImages(fname,mask)

% fname: output file of pvs_quantify.m
% m: mask for regions:  ind_bc (pvs that intersect centrum ovale and subcortical); ind_b (pvs's for subcortical); and ind_c (pvs's for centrum ovale)


    
    a=load(fname);
    
    d=double(a.c>0);
    if exist('mask','var')
    m=load(mask);
    
  %  m=[m.ind_bc,m.ind_b];
    
  %  for k=1:length(m)
   %   d(a.c==m(k))=0;
   % end
    
   
   %%{ 
    for k=1:length(m.ind_bc)
      d(a.c==m.ind_bc(k))=3;
    end
    
    for k=1:length(m.ind_b)
      d(a.c==m.ind_b(k))=2;
    end
    
    d=(d==1);
    
    end
        
    %}
    
    
   d=clusterize2(d,4);
   d2=zeros(size(d));
   d2b=d2;
%%{
for j=2%:length(a.path)
  for k=1:length(a.ind{j})
      
     pos=ind2subb(size(d),a.ind{j}(k));
     
     if d(pos(1),pos(2),pos(3))>0
       %  tmp=a.vclust{j}(k)/length(a.path{j})*254;
         tmp=a.vclust{j}(k)/20*255;
         if tmp>255
             tmp=255;
         end
         d2(pos(1),pos(2),pos(3))=tmp;
         d2b(pos(1),pos(2),pos(3))=a.vclust{j}(k);
         
     end
      
  end
  
  for k=1:length(a.path{j})
      
     pos=ind2subb(size(d),a.ind{j}(a.path{j}(k)));
     
     if d(pos(1),pos(2),pos(3))>0
         d2(pos(1),pos(2),pos(3))=256;
     end
      
  end
    
end
d2_3d=d2;
%}
%d=d(67:end-48,25:end-29,1:end-25); %testz
d2=d2(67:end-48,25:end-29,1:end); %HVB237, 238, 236, 239, 240

figure;
is=size(d2,3);
ap=size(d2,1);
lr=size(d2,2);
m=zeros(ap+is,lr+ap);

m(1:ap,1:lr)=Opac_mask(d2,3);
m(ap+1:ap+is,1:lr)=flipud(Opac_mask(d2,1)');
m(ap+1:ap+is,lr+1:lr+ap)=flipud(Opac_mask(d2,2)');
%{
m(1:ap,1:lr)=sum(d,3);
m(ap+1:ap+is,1:lr)=flipud(squeeze(sum(d,1))');
m(ap+1:ap+is,lr+1:lr+ap)=flipud(squeeze(sum(d,2))');
%}
cm=[0,0,0;jet(254);1,1,1];
writeanalyze(d2_3d,'pvs_path2');

writeanalyze(d2b,'pvs_path2b');

imshow(m,cm);
set(gcf,'Position',[51   191   977   730]);
%savetiffc(['PVS3D_',fname{i}]);



