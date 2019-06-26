function res = pvs_age_dependence(fname_list,age,thr,sym)
% res = pvs_stat(fname_list,thr,plot_sym)
%
% thr: the minimum volume for a pvs.

dname=fileparts(fileparts(fname_list{1}));

savename=fullfile(dname,sprintf('results_pvs_ageDependence_%dsub.mat',length(fname_list)));

if ~exist(savename,'file')
    
    
    diam_mean=zeros(length(fname_list),1);
    volume_mean=zeros(length(fname_list),1);
    length_mean=zeros(length(fname_list),1);
    
    cdiam=cell(length(fname_list),1);
    cvolume=cell(length(fname_list),1);
    clength=cell(length(fname_list),1);
    
    cdiam_pervessel=cell(length(fname_list),1);
    
    
    pvs_num = zeros(length(fname_list),1);
    
    
    for i=1:length(fname_list)
        
        
        a=load(fname_list{i});

        itmp=setdiff(1:length(a.l),find(a.v<thr));
        
        d=[];
        d_pervessel=[];
        
        for j2=1:length(itmp)
            d=[d(:);a.dpath{itmp(j2)}(:)];
            d_pervessel(j2)=mean(a.dpath{itmp(j2)}(:));
        end
        
        diam_mean(i)=mean(d);
        volume_mean(i)=mean(a.v);
        length_mean(i)=mean(a.l);
        pvs_num(i)=length(a.l);
        
        cdiam{i}=d;
        cvolume{i}=a.v(itmp);
        clength{i}=a.l(itmp);
        cdiam_pervessel{i}=d_pervessel;
    end
    
    fprintf('Mean diameters in subjects :\n');
    disp(diam_mean);
    fprintf('Group mean: \n');
    disp(mean(diam_mean,1));
    fprintf('STD: \n');
    disp(std(diam_mean,[],1));
    
    %%
    save(savename,'age','diam_mean');
    
else
    
    
    load(savename);
end

    hold on;
    myScatter(age,diam_mean,'Age (yrs)','Mean Diameter (mm)',2,sym);
    
    