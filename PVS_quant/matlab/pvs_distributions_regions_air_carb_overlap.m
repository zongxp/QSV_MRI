function [x2,hist_air,raw_air,hist_carb,raw_carb]=pvs_distributions_regions_air_carb_overlap(fname_list_air,fname_list_carb,ciPVS_match)

% input:
% file should contain:
% dpath, ind_n (include ambiguous PVS), ind_n2 (exclude ambiguous PVS), l, v
% ind_n and ind_n2: 1*nregions cells;
% ciPVS_match: PVS indices to include in the ouput; 1*nsub cells;

% output:
% x2: 1*3 cell: bin positions
% hist_feature: nregion*3 cells; contains counts for the histogram;
% nsub*nbin
% raw_feature:  nregion*3*nsub cells: raw values; not converted to histograms


ncase=length(fname_list_air);
nregions=get_nregion(fname_list_air);

hist_feature=cell(nregions,3,2);  % nregion regions; 3 features; histograms; 2 gas conditions    
raw_feature = cell(nregions,3,2);
use_amb=true; % when ambiguous, exclude or not.

for i=1:ncase
    %m=load(sprintf('PVS_region_check_%s_4lobes.mat',regioncheck_list{i}));
    disp(fname_list_air{i});
    disp(fname_list_carb{i});
    
    if isempty(fname_list_air{i}) || ~exist(fname_list_air{i},'file') ||isempty(fname_list_carb{i}) || ~exist(fname_list_carb{i},'file')
        continue;
    end
    
    a=load(fname_list_air{i});   
    a2=load(fname_list_carb{i});
    
    
    lbin=[0:2:80,inf]*0.4;  %length
    vbin=[0:4:300,inf]*0.4^3; % volume
    dbin=[0:0.2:10,inf]*0.4;
        
    x2{1}=(lbin(1:end-1)+(lbin(2)-lbin(1))/2);
    x2{2}=(vbin(1:end-1)+(vbin(2)-vbin(1))/2);
    x2{3}=(dbin(1:end-1)+(dbin(2)-dbin(1))/2);
    
    for j=1:nregions%length(m.ind_n2)
        

         if use_amb
            [i_air,i_carb,i_air2,i_carb2]=get_ind(a.l,a2.l,a.ind_n{j},a2.ind_n{j},ciPVS_match{i});
         else
            [i_air,i_carb,i_air2,i_carb2]=get_ind(a.l,a2.l,a.ind_n2{j},a2.ind_n2{j},ciPVS_match{i});
         end

        [hist_air{j,1}(i,:),hist_air{j,2}(i,:),raw_air{j,1,i},raw_air{j,2,i}]=get_hist_raw(a,i_air,lbin,vbin,i_air2);
        [hist_carb{j,1}(i,:),hist_carb{j,2}(i,:),raw_carb{j,1,i},raw_carb{j,2,i}]=get_hist_raw(a2,i_carb,lbin,vbin,i_carb2);
             
        [hist_air{j,3}(i,:),raw_air{j,3,i}]=get_hist_raw_diam(a,i_air,dbin);
        [hist_carb{j,3}(i,:),raw_carb{j,3,i}]=get_hist_raw_diam(a2,i_carb,dbin);
        
    end
    
end

function [hist_l,hist_v,raw_l,raw_v]=get_hist_raw(a,i,lbin,vbin,i2)

if isempty(i)
    hist_l=0*lbin;
    hist_v=0*vbin;
else
    hist_l=histc(a.l(i),lbin);
    hist_v=histc(a.v(i),vbin);
end
raw_l=a.l(i2);  % include PVS with length < 0.8 mm;
raw_v=a.v(i2);

function [hist_d,raw_d]=get_hist_raw_diam(a,i,dbin)

        d=[];
        for j2=1:length(i)
            d=[d(:);a.dpath{i(j2)}(2:end-1)'];
        end
        if ~isempty(d)
            hist_d=histc(d,dbin);
        else
            hist_d=0*dbin;
        end
        raw_d=d';  %do not include PVS with length<0.8 mm

function [i_air,i_carb,i_air2,i_carb2]=get_ind(l_air,l_carb,ind_air,ind_carb,ind_match)
% l_air, l_carb: length of PVS; 1*n 
% ind_air, ind_carb: indices of PVS for the roi; 1*n2 where n2<=n
% ind_match: matched indices; n3*2 where n3<=n 
% output: i_air, i_carb: only include data for length >= 0.8
%         i_air2, i_carb2: include all data.

i_air=[];
i_carb=[];
i_air2=[];
i_carb2=[];

for i=1:size(ind_match,1)

    if l_air(ind_match(i,1))>=0.8 && l_carb(ind_match(i,2))>=0.8 && any(ind_match(i,1)==ind_air) && any(ind_match(i,2)==ind_carb)  
      i_air(end+1)=ind_match(i,1);
      i_carb(end+1)=ind_match(i,2);
    end
    
    if  any(ind_match(i,1)==ind_air) && any(ind_match(i,2)==ind_carb)   
      i_air2(end+1)=ind_match(i,1);
      i_carb2(end+1)=ind_match(i,2);
    end
    
end


function nregions=get_nregion(fname_list)


ncase=length(fname_list);
for i=1:ncase
   %m=load(sprintf('PVS_region_check_%s_4lobes.mat',regioncheck_list{i}));
%   disp(fname_list{i});
   if isempty(fname_list{i}) || ~exist(fname_list{i},'file')
       continue;
   end
ind_n=ri(fname_list{end},'','','ind_n');
nregions=length(ind_n);
break;
end


