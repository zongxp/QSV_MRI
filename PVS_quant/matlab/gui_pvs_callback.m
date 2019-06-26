function gui_pvs_callback(params,sbutton)

if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end

    files=get_fpattern(params,'mask files');
    files_pattern=get(params,'mask files');
 
    T2file_pattern=get(params,'T2 files');
    
    prob_files=get_fpattern(params,'prob files',1);
    files=str2cell(files);

    
    prob_thr=get(params,'prob thr');
    clust_size=get(params,'clust size thr'); 
    
  
if strcmp(sbutton,'process')

    for i=1:length(files)
        if ~isempty(T2file_pattern)
                      
         T2=wildcard2id(filename(files_pattern),filename(files{i}),filename(T2file_pattern));
               
         
         fprintf('process_pvs(''%s'', ''%s'') \n',filename(files{i}),fullfile(fileparts(T2file_pattern),T2));
      
         process_pvs(files{i},fullfile(fileparts(T2file_pattern),T2));
        else
         fprintf('process(''%s'') \n',files{i});   
         process_pvs(files{i});
        end
      
    end
%     gui_pvs_callback(params,'show results');
elseif strcmp(sbutton,'gen mask')
    
   
    for i=1:length(prob_files)
        
        
        nii_name=strtok2(prob_files{i},'.');
        prefix=strtok2(nii_name,'.');
        prefix=filename_prefix(prefix,'mask_');
        
        if ~exist([prefix,'.nii.gz'],'file')
            d=load_untouch_niigz(prob_files{i});
            
            [m,nvox_sort]=clusterize2(d.img>prob_thr,clust_size);
            
            d.img=m>0;
            fprintf('%s - Number of PVSs found = %d\n',filename(prefix),max(m(:)));
            save_untouch_niigz(d,prefix,true);
        end
    end
    
elseif strcmp(sbutton,'show results')  
  
%     
% elseif strcmp(sbutton,'gen movie')  
%     fT2=get(params,'T2 image');
%    fmask=get(params,'mask files');
%     im=ri(fT2);
%     intmax=get(params,'max intensity');
%     m=ri(fmask);
%     
%     delayS = get(params,'delay (s)');
%    im2=scale2n(im,100,[0,intmax]);
%   % im2=repmat2(im2,2);
%   % m=repmat2(m,2);
%   % m=bwmorph3d(m,'remove');
%    [im3,cm]=combine_over_under(im2,m,gray(100),[1,0,0],m>0);
%    
%    im4=cat(4,im2,im3);
%    im4=permute(im4,[1,2,4,3]);
%    im4=reshape(im4,[size(im4,1),size(im4,2),2*size(im4,4)]);
%    vol2gif(im4,cm,strtok(fT2,'.'),delayS)
   
elseif strcmp(sbutton,'gen bsub')
    
      files=get_fpattern(params,'mask files');
    
      files=str2cell(files);
      fid=fopen('bsub_pvs.sh','w');
    for i=1:length(files)
      
        fname=filename(files{i});
        fprintf(fid,'bsub -M 8 -o logpvs_%s.%%J matbgk "process_pvs(''%s'',%3.2f)" logpvs_%s\n',strtok(fname,'.'),fname,vox_size,strtok(fname,'.'));
        
    end
    fclose(fid);
 
    
elseif strcmp(sbutton,'pvs mip')
    fmask=get(params,'mask files');
    m=ri(fmask);
    m=m>0;
   
  %  m=Skeleton3D(m);
    alpha=0.9;
    
    dim=get(params,'dimension');
    
    m=adjust_dim(m,dim);
   
    m2=cumsum(m,3);
    m3=alpha.^(m2-1).*m;
    m4=squeeze(sum(m3,3)); 
    figure;imshow(m4,[]);
    
 %   savetiffc(prefix);


elseif strcmp(sbutton,'pvs rotate gif')
    fmask=get(params,'mask files');
    m=ri(fmask);
    m=m>0;
   
    dim=get(params,'dimension');

     m=adjust_dim(m,dim);
    
    prefix=filename_prefix(strtok2(fmask,'.'),'rot_');
    vol2gif_rotate(m,[1,3],30,0.2,prefix); 
    
 %   savetiffc(prefix);

end

disp('gui_pvs_callback done');


function m=adjust_dim(m,dim)

for i=1:length(dim)
    if dim(i)<0
        m=flip(m,abs(dim(i)));
    end
end
m=permute(m,abs(dim));  %make the first dimension AP, second LR, third, IS


function res=wildcard2id(wildcard,matched,wildcard2)
i=find(wildcard=='*');
if length(i)>1
    error('not support more than 1 *');
elseif isempty(i)
    
    res=wildcard2;
    return;
end

matched_pattern=wildcard_matched_pattern(wildcard,matched);

i2=find(wildcard2=='*');
if length(i2)>1
    error('not support more than 1 *');
end

res=[wildcard2(1:i2-1),matched_pattern,wildcard2(i2+1:end)];



function res=wildcard_matched_pattern(wildcard,matched)
i=find(wildcard=='*');

if length(i)>1
    error('not support more than 1 *');
end
n=length(wildcard)-i;
res=matched(i:end-n);




