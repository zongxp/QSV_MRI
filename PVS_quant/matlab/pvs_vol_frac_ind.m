function pvs_vol_frac_ind(fname,t2)

%fname: should be the *_Dv_stat.mat file
%t2: T2 images

fprintf('run pvs_vol_frac_ind(%s,%s)\n',fname,filename(t2));

if strcmp(fname(end-3:end),'.mat')
    fname=fname(1:end-4);
end

t2=ri_d1(t2);
a=load(fname);
voxsize=a.voxsize;
ind=a.ind;
c=a.c;

vf_ind=zeros(1,length(ind));


c_norm=zeros(1,length(ind)); % normalized contrast relative to WM

nvox=zeros(1,length(ind));

for i=1:length(ind)

 ind_ext=m_ext_calc(size(c),ind{i});
 
 cext=c(ind_ext);
 ind_ext(cext>0)=[];

  s_wm=mean(t2(ind_ext));
  
  nvox(i)=length(ind{i});
  spvs=mean(t2(ind{i}));
  
  vf_ind(i)=(spvs/s_wm-1)*nvox(i)*prod(voxsize); % =v_pvs*(s_pvs/s_wm-1);
  c_norm(i)=(spvs/s_wm-1);
  
    
end

prefix=fname;
save([prefix,'_vfind.mat'],'vf_ind','voxsize','nvox','c_norm');
 copy_orient(fname,[prefix,'_vfind.mat']);
return;
%%

function ind2=m_ext_calc(sz,ind)


l=2;
ind2=zeros(length(ind)*(2*l+1)^3,1);
count=0;
for i=1:length(ind)

    ijk=ind2subb(sz,ind(i));

                if ijk(1)-l<1 || ijk(1)+l>sz(1) || ijk(2)-l<1 || ijk(2)+l>sz(2) ||ijk(3)-l<1 || ijk(3)+l>sz(3)    
                    continue;
                end
                
    for i1=-l:l
        for i2=-l:l
            for i3=-l:l
                
                
                count=count+1;
               ind2(count)=ind(i)+i1+i2*sz(1)+i3*(sz(1)*sz(2));
                
            end
        end
    end
    
end
ind2=setdiff(ind2,ind);

ind2(ind2==0)=[];





