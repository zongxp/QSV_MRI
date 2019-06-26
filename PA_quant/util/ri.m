function [d,dim]=ri(fname,ns,do_sort,var_name)
%ri(fname,ns,do_sort,field_name)
% ns: number of slices in the image; when empty; will be set to the
% maximum possible.
% var_name: used for reading .mat file; which variable to read, default:roi

if ~isa(fname,'char')  % it is a matrix;
    d=fname;
    dim=ones(1,ndims(d));
    return;
end

if ~exist('do_sort','var')
    do_sort=true;
end

[fname2,r] = strtok(fname,'[');
sel=[];
if ~isempty(r)
    r =strtok(r(2:end),']');
    [r,r2] = strtok(r,'.');
    if isempty(r2)
        sel= str2double(r)+1;
    else
        sel = str2double(r)+1:str2double(r2(3:end))+1;
    end
end
dim=0;



[dname,fname3,ext]=fileparts(fname2);

if strcmp(ext,'.HEAD') || exist([fname2,'.HEAD'],'file')
    
    d=BrikLoadf(fname);
elseif isdir(fname2) || strcmp(ext,'.dcm') || strcmp(ext,'.IMA') || isempty(ext)
    
    if exist('ns','var')
        d=dcm2mat(dname,fname3,ext,do_sort,ns);
    else
        d=dcm2mat(dname,fname3,ext,do_sort);
    end
elseif strcmp(ext,'.nii')
    d=load_untouch_nii(fname2);
    d=d.img;
    
elseif strcmp(ext,'.gz')
    
    d=load_untouch_niigz(fname2);
    d=d.img;
    
    
elseif strcmp(ext,'.hdr') || exist([fname2,'.hdr'],'file')
    [d,dim]=readanalyze(fname);
elseif strcmp(ext,'.xml') || exist([fname2,'.xlm'],'file')
    d=read_data_amide(fname3(10:end),dname);
    if ndims(d)==5
        sz=size(d);
        d=reshape(d,sz([1,2,3,5]));
    end
elseif strcmp(ext,'.mat') || exist([fname2,'.mat'],'file')
    tmp=whos('-file',fname);
    
    for i=1:length(tmp)
        nm{i}=tmp(i).name;
    end
    
    if exist('var_name','var')
        
        if any(strcmp(var_name,nm))
            tmp=load(fname,var_name);
            d=eval(sprintf('tmp.%s',var_name));
        else
            error('Field name %s not exist',var_name);
        end
    else
        if length(nm)==1
            tmp=load(fname);
            d=eval(sprintf('tmp.%s',nm{1}));
            
            if isstruct(d) % generated in load_untouch_niigz;
                d=d.img;
            end
            
        else
            
            
            
            fprintf('Field names in file: %s.\n',fname);
            disp(nm);
            fn=inputdlg({'select the field name:' },'',1,{sprintf('%s',nm{1})});
            
            tmp=load(fname,fn{1});
            d=eval(sprintf('tmp.%s',fn{1}));
        end
        
        
    end
    
    dim=[1,1,1];
else
    fname=strtok(fname,'.');
    d=rdSdt(fname);
    
    d=flipdim(d,1);
end


if nargout==0
    d=[];
end

function d=dcm2mat(dname,fname2,ext,do_sort,ns)

warning('off','all');

if isdir(fname2)
    dname=fname2;
elseif isempty(dname)
    dname=pwd;
end

dir_str=dir(fullfile(dname,['*',ext]));

d=[];
spos=[];
i2=0;
SliceLocation=[];
% [voxsize,center,orient]=dcmDimCenter2(dname);
for i=1:length(dir_str)  %skip '.' and '..'
    
    
    if mod(i,10)==0
        disp(i);
    end
    
    if isdir(fullfile(dname,dir_str(i).name))
        continue;
    end
    
    
    %  disp([i,size(tmp)]);
    try
        h=dicominfo(fullfile(dname,dir_str(i).name));
    catch
        continue;
    end
    i2=i2+1;
    tmp=dicomread(fullfile(dname,dir_str(i).name));
    spos(:,i2)=h.ImagePositionPatient;
    SliceLocation(i2)=h.SliceLocation;
    %  disp([i,size(tmp)]);
    if ~isempty(d)&& size(tmp,2)==size(d,1) && size(tmp,1)==size(d,2) && size(d,1)~=size(d,2)
        d(:,:,i2)=tmp';   %for some dicom images, the dimension switches.
    else
        d(:,:,i2)=tmp;
    end
    
    
end

if ~exist('ns','var')
    
    ns=size(d,3);
    %%{
    if size(d,3)>1
        ns=inputdlg({sprintf('Number of slices (max=%d)',size(d,3))},'',1,{sprintf('%d',size(d,3))});
        ns=str2num(ns{1});
    else
        ns=1;
    end
    
    %}
end

if isempty(ns)
    ns=size(d,3);
end

if mod(size(d,3),ns)~=0
    error('number of slices not possible');
end
d=reshape(d,[size(d,1),size(d,2),ns,size(d,3)/ns]);

if do_sort && ns>1
    spos=reshape(spos,[3,ns,size(d,4)]);
    [postmp,ind]=sort(unique(SliceLocation));
    d=d(:,:,ind,:);
    spos=spos(:,ind,1);
end

if ns==1
    spos=spos(:,1);
end

[voxsize,center,orient,pos,rotmat]=dcmDimCenter(spos,h.ImageOrientationPatient,h.PixelSpacing,size(d,2),size(d,1));

if size(spos,2)==1
    voxsize(3)=readdPar(dname,'SliceThickness');
end

save([dname,'.mat'],'d','voxsize','center','orient','pos','rotmat');
 warning('on','all');
