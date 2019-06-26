function PVS_orient(fname)

if strcmp(fname(end-3:end),'.mat')
    fname=fname(1:end-4);
end

a=load(fname);
nvox_path=a.nvox_path;
ind=a.ind;
sub=a.sub;

c=a.c;
npath=length(nvox_path);
lthr=0.8; %lthreshold = 3 mm
sz=size(c);
norm=zeros(3,npath);  % normalized direction of the whole pathes
norm_seg=cell(1,npath); %average direction of 5 neighboring voxels in the paths

cnorm=zeros([3,sz]);
cnorm_seg=NaN([3,sz]);


for i=1:npath
 
    posc=sub{i};  
    posc=a.rotmat*posc'.*a.voxsize(:);
    
    orient=pca(posc');
    
    norm(:,i)=orient(:,1);
    
  if 0  
    cnorm(:,ind{i})=repmat(norm(:,i),[1,length(ind{i})]);
  end
    pos=posc(:,1:nvox_path(i));

    norm_seg{i}=NaN(3,size(pos,2));
    
    for j=1:size(pos,2)
        
        if size(pos,2)>=5
            if j-2>0 && j+2<size(pos,2)
                orient=pca(pos(:,j-2:j+2)');
            elseif j<=2
                orient=pca(pos(:,1:5)');
            else
                orient=pca(pos(:,end-4:end)');
            end
                
               norm_seg{i}(:,j)=orient(:,1);
                 
        end
        
    end

    if 0 %no need to save
      for j=1:length(ind{i})
            
          dist=sos(posc(:,j)-pos,1);
          
          [mdist,ind_min]=min(dist);
         try
           cnorm_seg(:,ind{i}(j)) = norm_seg{i}(:,ind_min);
         
         catch
            disp(''); 
         end
         % end
      end
      
    
    end
      
end
%cnorm=permute(cnorm,[2,3,4,1]);
%cnorm_seg=permute(cnorm_seg,[2,3,4,1]);

orient=get_orient(fname);
%orient.cnorm=cnorm;
%orient.cnorm_seg=cnorm_seg;
orient.norm=norm;
orient.norm_seg=norm_seg;
save([fname,'_orient'],'-struct','orient');





