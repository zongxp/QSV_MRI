function [x2,hist_feature,raw_feature,roi_nvox,voxsize]=pvs_distributions_regions(fname_list,fname_list_vfind,cPVS_ind)

% input:
% file should contain:
% dpath, ind_n (include ambiguous PVS), ind_n2 (exclude ambiguous PVS), l, v
% ind_n and ind_n2: 1*nregions cells;
% cPVS_ind: PVS indices to include in the ouput; 1*nsub cells;

% output:
% x2: 1*3 cell: bin positions
% hist_feature: nregion*3 cells; contains counts for the histogram;
% nsub*nbin
% raw_feature:  nregion*5*nsub cells: raw values; not converted to
% histograms; 4 mean diameter; 5 max diameter;

ncase=length(fname_list);

nregions=get_nregion(fname_list);

hist_feature=cell(nregions,3);  % nregion regions; 3 features; histograms    
raw_feature = cell(nregions,3);
use_amb=true; % when ambiguous, exclude or not.

roi_nvox=zeros(ncase,nregions);
for i=1:ncase
    %m=load(sprintf('PVS_region_check_%s_4lobes.mat',regioncheck_list{i}));
    disp(fname_list{i});
    if isempty(fname_list{i}) || ~exist(fname_list{i},'file')
        continue;
    end
    a=load(fname_list{i});
    
    if ~isempty(fname_list_vfind)
     a_vf=load(fname_list_vfind{i});
    end
    
    voxsize=a.voxsize;
    
    roi_nvox(i,:)=a.roi_nvox;
       
    lbin=[0:2:10,14:4:80,inf]*a.voxsize(1);  %length
    vbin=[0:3:15,20:5:300,inf]*prod(a.voxsize); % volume
    dbin=[0:0.2:5,inf]*a.voxsize(1);
        
    %x2{1}=(lbin(1:end-1)+lbin(2:end))/2;
    %x2{2}=(vbin(1:end-1)+vbin(2:end))/2;
    %x2{3}=(dbin(1:end-1)+dbin(2:end))/2;
    
    x2{1}=lbin;
    x2{2}=vbin;
    x2{3}=dbin;
    
    
    for j=1:nregions%length(m.ind_n2)
        
        if use_amb
            itmp=a.ind_n{j};
        else
            itmp=a.ind_n2{j};
        end
        exc=[];
        exc2=[];
        for j2=1:length(itmp)
            if a.l(itmp(j2))<0.8 || (exist('cPVS_ind','var') && ~any(itmp(j2)==cPVS_ind{i})) %|| a.v(itmp(j2))*vox_size^3>25
                exc(end+1)=j2;
            end
            if  exist('cPVS_ind','var') &&~any(itmp(j2)==cPVS_ind{i}) %|| a.v(itmp(j2))*vox_size^3>25
                exc2(end+1)=j2;
            end
        end
        
        itmp2=itmp;
        itmp(exc)=[];
        itmp2(exc2)=[];
        
        
        hist_feature{j,1}(i,:)=histc2(a.l(itmp),lbin);
        
        hist_feature{j,2}(i,:)=histc2(a.v(itmp),vbin);
        raw_feature{j,1,i}=a.l(itmp2);  % include PVS with length < 0.8 mm;
        raw_feature{j,2,i}=a.v(itmp2);
        
        d=[];
        dmean=[];
        dmax=[];
        for j2=1:length(itmp)
            d=[d(:);a.dpath{itmp(j2)}(2:end-1)'];
            dmean(j2)=mean(a.dpath{itmp(j2)}(2:end-1));
            
            if isempty(a.dpath{itmp(j2)}(2:end-1))
             dmax(j2)=NaN;
            else
             dmax(j2)=max(a.dpath{itmp(j2)}(2:end-1));
            end
            
        end
%%        
if ~isempty(fname_list_vfind)
    vf_ind=[];
    for j2=1:length(itmp)
        vf_ind(j2)=a_vf.vf_ind(itmp(j2));
    end
    raw_feature{j,6,i}=vf_ind;
end
%%        
        hist_feature{j,3}(i,:)=histc2(dmean,dbin);       
        raw_feature{j,3,i}=d';  %do not include PVS with length<0.8 mm
        raw_feature{j,4,i}=dmean;
        raw_feature{j,5,i}=dmax;
      
        
        
    end
    
end


function nregions=get_nregion(fname_list)


ncase=length(fname_list);
for i=1:ncase
   %m=load(sprintf('PVS_region_check_%s_4lobes.mat',regioncheck_list{i}));
%   disp(fname_list{i});
   if isempty(fname_list{i})
       continue;
   end
ind_n=ri(fname_list{end},'','','ind_n');
nregions=length(ind_n);
break;
end


