function res = pvs_stat(fname_list,thr,plot_sym)
% res = pvs_stat(fname_list,thr,plot_sym)
%
% thr: the minimum volume for a pvs.


dname=fileparts(fileparts(fname_list{1}));

savename=fullfile(dname,sprintf('results_pvs_stat_%dsub.mat',length(fname_list)));

if ~exist(savename,'file')
    
    if ~exist('plot_sym','var')
        plot_sym='bo';
    end
    
    gap=0:0.05:1;
    
    xmn=zeros(length(fname_list),size(gap,2)-1);  %scaled position
    ymn=zeros(length(fname_list),size(gap,2)-1);  % diameter
    
    
    prct=[25,50,75];
    res_iqr=zeros(length(fname_list),length(prct),1,3);  %  1 regions.3 features;
    diam_mean=zeros(length(fname_list),1);
    volume_mean=zeros(length(fname_list),1);
    length_mean=zeros(length(fname_list),1);
    
    cdiam=cell(length(fname_list),1);
    cvolume=cell(length(fname_list),1);
    clength=cell(length(fname_list),1);
    
    cdiam_pervessel=cell(length(fname_list),1);
    
    
    pvs_num = zeros(length(fname_list),1);
    
    y2=cell(1,3);  % length; volume; diameter distributions
    
    for i=1:length(fname_list)
        
        disp(fname_list{i});
        a=load(fname_list{i});
        
        
        xall=[];
        yall=[];
        
        
        for j=1:length(a.vpath)
            
            x=a.lnorm{j}; %linspace(0,1,length(a.vpath{j}));
            y=a.dpath{j};
            
            xall=[xall,x];
            yall=[yall,y];
            
        end
        
        for j=1:size(gap,2)-1
            if j==1
                ind=xall>gap(j)&xall<gap(j+1);
            else
                ind=xall>=gap(j)&xall<gap(j+1);
            end
            xmn(i,j)=mean(xall(ind));
            ymn(i,j)=mean(yall(ind));
            
        end
        
        %%{
        
        rbin=[0:0.1:2,inf]; % radius
        abin=[0:0.3:10,inf]; % cross-sectional area
        
        
        exc=a.l<=2;
        
        
        lbin=[0:2:80,inf]*0.4;  %length
        vbin=[0:4:300,inf]*0.4^3; % volume
        dbin=[0:0.2:10,inf]*0.4;
        
        
        x2{1}=(lbin(1:end-1)+(lbin(2)-lbin(1))/2);
        x2{2}=(vbin(1:end-1)+(vbin(2)-vbin(1))/2);
        x2{3}=(dbin(1:end-1)+(dbin(2)-dbin(1))/2);
        
        
        
        itmp=setdiff(1:length(a.l),find(a.v<thr));
        y2{1}(:,i)=histc(a.l(itmp),lbin);
        y2{2}(:,i)=histc(a.v(itmp),vbin);
        res_iqr(i,:,1,1)=prctile(a.l(itmp),prct);
        
        res_iqr(i,:,1,2)=prctile(a.v(itmp),prct);
        
        
        
        d=[];
        d_pervessel=[];
        
        for j2=1:length(itmp)
            d=[d(:);a.dpath{itmp(j2)}(:)];
            d_pervessel(j2)=mean(a.dpath{itmp(j2)}(:));
        end
        y2{3}(:,i)=histc(d,dbin);
        res_iqr(i,:,1,3)=prctile(d,prct);
        
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
    save(savename);
    
else  
    load(savename);
end   
    subplot(2,2,4);
    hold on;
    
    errorbar(mean(xmn,1),mean(ymn,1),std(ymn,[],1)/sqrt(length(fname_list)),plot_sym,'LineWidth',1.5);
    box on;
    
    set(gca,'FontSize',14);
    xlabel('Normalized Distance');
    ylabel('Diameter (mm)');
    ylim([0.58 0.9]);
    set(gca,'XTick',0:0.2:1);
    %ylm=[0,0.3;0,0.2;0,0.2;0,0.3];
    %xlm=[0,30;0,15;0,0.6;0,20];
    ylm=[0,0.22;0,0.2;0,0.2];
    xlm=[0,20;0,8;0,2];
    
    mres_iqr=(mean(res_iqr,1));
    eres_iqr=(std(res_iqr,[],1));
    xlb={'Length (mm)','Volume (mm^3)','Diameter (mm)'};
    for i=1:3
        
        
        disp(xlb{i});
        fprintf('Interquitile range: ');
        fprintf('%10d',prct);
        fprintf('\n');
        fprintf('WM                : ');
        fprintf('%10.2e',mres_iqr(1,:,1,i));
        fprintf('\n');
        fprintf('                    ');
        fprintf('%10.2e',eres_iqr(1,:,1,i));
        fprintf('\n');
        
        subplot(2,2,i);
        %errorbar(x{i},mean(y{i}(1:end-1,:),2),std(y{i}(1:end-1,:),[],2),'ko-','LineWidth',1.5);
        %hold on;
        hold on;
        %disp(xlb{i});
        
        y=y2{i};
        
        yn=sum(y,1);
        yn=repmat(yn,[size(y,1),1]);
        y=y./yn;
        [tmp,ind_max]=max(y,[],1);
      
        errorbar(x2{i},mean(y(1:end-1,:),2),std(y(1:end-1,:),[],2)/sqrt(length(fname_list)),plot_sym,'LineWidth',1.5);
      
        if i==2
            set(gca,'Xtick',0:2:10);
        elseif i==3
            set(gca,'Xtick',0:0.4:2);
            set(gca,'XMinorTick','On');
        end
        
        set(gca,'FontSize',14);
        xlabel(xlb{i});
        ylabel('Probability');
        %xlim([0,x{i}(end-2)]);
        ylim(ylm(i,:));
        xlim(xlm(i,:));
        box on;
        
    end
    
%savetiffc(sprintf('pvs_stat_%dsub',length(fname_list)));

