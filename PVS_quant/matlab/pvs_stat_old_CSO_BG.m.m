function pvs_stat(fname_list,regioncheck_list,vox_size)


%dname='.';
n=[];
n_bc=[];
n_b=[];

    gap=0:0.05:1;
    
    xmn=zeros(length(fname_list),size(gap,2)-1);
    ymn=zeros(length(fname_list),size(gap,2)-1);
   
    xmn_th=zeros(length(fname_list),size(gap,2)-1,2);
    ymn_th=zeros(length(fname_list),size(gap,2)-1,2);
    
    prct=[25,50,75];
    res_iqr=zeros(length(fname_list),length(prct),2,3);  %  2 regions.3 features;
    diam_mean=zeros(6,4);  %% 6 subjects; 3 regions (wm,sc, intersect both, whole brain);
    
for i=1:length(fname_list)
   m=load(sprintf('PVS_region_check_%s.mat',regioncheck_list{i}));

    a=load(sprintf('PVSLength_%s_Curv_Dv_stat',fname_list{i}));
    

xall=[];
yall=[];

xall_sc=[];
yall_sc=[];

yall_th=cell(1,2);
xall_th=cell(1,2);

    for j=1:length(a.vpath)
       
          
       if m.path_dir2(j)==1     
           itmp=m.ind_n2{1};
           if ~any(itmp==j)
               continue;
           end
        x=a.lnorm{j}; %linspace(0,1,length(a.vpath{j}));
        y=a.dpath{j};
        
        xall=[xall,x];
        yall=[yall,y];
       elseif m.path_dir2(j)==-1
           itmp=m.ind_n2{1};
           if ~any(itmp==j)
               continue;
           end
           
        x=1-a.lnorm{j};
        y=a.dpath{j};
        
        xall=[xall,x];
        yall=[yall,y];
       else
           
          
          itmp=m.ind_n2{2};
           if ~any(itmp==j)
               continue;
           end
          
           xall_sc=[xall_sc,a.lnorm{j}];
           yall_sc=[yall_sc,a.dpath{j}];
           
           continue;
       end
       
 
        if m.path_th(j)>pi/4
            yall_th{2}=[yall_th{2},y];
            xall_th{2}=[xall_th{2},x];
        else
            yall_th{1}=[yall_th{1},y];
            xall_th{1}=[xall_th{1},x];
        end
            
    end
    
    for j=1:size(gap,2)-1
       if j==1
           ind=xall>gap(j)&xall<gap(j+1);
       else
        ind=xall>=gap(j)&xall<gap(j+1);
       end
       xmn(i,j)=mean(xall(ind));
       ymn(i,j)=mean(yall(ind));  
       
       if j==1
        ind=xall_sc>gap(j)&xall_sc<gap(j+1);
       else
        ind=xall_sc>=gap(j)&xall_sc<gap(j+1);
       end
       
       xmn_sc(i,j)=mean(xall_sc(ind));
       ymn_sc(i,j)=mean(yall_sc(ind)); 
       
      for k=1:2 
       if j==1
        ind=xall_th{k}>gap(j)&xall_th{k}<gap(j+1);
       else
        ind=xall_th{k}>=gap(j)&xall_th{k}<gap(j+1);
       end
       xmn_th(i,j,k)=mean(xall_th{k}(ind));
       ymn_th(i,j,k)=mean(yall_th{k}(ind)); 
      end
       
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
d_brain=[];

for j=1:length(m.ind_n2)
    
    itmp=m.ind_n2{j};
    
    if ~isempty(itmp)
     y2{1,j}(:,i)=histc(a.l(itmp)*vox_size,lbin);
     y2{2,j}(:,i)=histc(a.v(itmp)*vox_size^3,vbin);
     res_iqr(i,:,j,1)=prctile(a.l(itmp)*vox_size,prct);
     
     res_iqr(i,:,j,2)=prctile(a.v(itmp)*vox_size^3,prct);
     
     
    else
      y2{1,j}(:,i)=zeros(length(lbin),1);
      y2{2,j}(:,i)=zeros(length(vbin),1);
    end
    

    d=[];
    for j2=1:length(itmp)        
        d=[d(:);a.dpath{itmp(j2)}(:)];
    end
    
    if ~isempty(d)
     y2{3,j}(:,i)=histc(d*vox_size,dbin);
     res_iqr(i,:,j,3)=prctile(d*vox_size,prct);
     
      
    else
     y2{3,j}(:,i)=zeros(length(dbin),1);
    end
    d_brain=[d_brain(:);d];
    diam_mean(i,j)=mean(d);
    
end

diam_mean(i,4)=mean(d_brain);
%}
end

fprintf('Mean diameters in subjects (sub*region (wm, sc, intersect both, whole brain)):\n');
disp(diam_mean*vox_size);
fprintf('Group mean: \n');
disp(mean(diam_mean*vox_size,1));
fprintf('STD: \n');
disp(std(diam_mean*vox_size,[],1));

%%
dmn=ymn*vox_size;


dmn_sc=ymn_sc*vox_size;


figure;
subplot(2,2,4);
hold on;
    
    errorbar(mean(xmn,1),mean(dmn,1),std(dmn,[],1)/sqrt(length(fname_list)),'or-','LineWidth',1.5);
    errorbar(mean(xmn_sc,1),mean(dmn_sc,1),std(dmn_sc,[],1)/sqrt(length(fname_list)),'bs-','LineWidth',1.5);
    
    set(gca,'FontSize',14);
    xlabel('Normalized Distance');
    ylabel('Diameter (mm)');
ylim([0.58 0.9]);
set(gca,'XTick',0:0.2:1);
%ylm=[0,0.3;0,0.2;0,0.2;0,0.3];
%xlm=[0,30;0,15;0,0.6;0,20];
ylm=[0,0.22;0,0.2;0,0.2];
xlm=[0,20;0,8;0,2];

xlb={'Length (mm)','Volume (mm^3)','Diameter (mm)'};
for i=1:3
    
    
mres_iqr=squeeze(mean(res_iqr,1));
eres_iqr=squeeze(std(res_iqr,[],1));
disp(xlb{i});
fprintf('Interquitile range: ');
fprintf('%10d',prct);
fprintf('\n');
fprintf('WM                : ');
fprintf('%10.2e',mres_iqr(:,1,i));
fprintf('\n');
fprintf('                    ');
fprintf('%10.2e',eres_iqr(:,1,i));
fprintf('\n');

fprintf('BG                : ');
fprintf('%10.2e',mres_iqr(:,2,i));
fprintf('\n');
fprintf('                    ');
fprintf('%10.2e',eres_iqr(:,2,i));
fprintf('\n');

subplot(2,2,i);
%errorbar(x{i},mean(y{i}(1:end-1,:),2),std(y{i}(1:end-1,:),[],2),'ko-','LineWidth',1.5);
%hold on;
 hold on;
sym={'ro-','bs-'};
%disp(xlb{i});
  for j=1:2
       y=y2{i,j};
if i==3 || i==1 || i==2

yn=sum(y,1);
yn=repmat(yn,[size(y,1),1]);
y=y./yn;
[tmp,ind_max]=max(y,[],1);






%fprintf('Peak position = %f (%f)\n',mean(x2{i}(ind_max)),std(x2{i}(ind_max)));


errorbar(x2{i},mean(y(1:end-1,:),2),std(y(1:end-1,:),[],2)/sqrt(length(fname_list)),sym{j},'LineWidth',1.5);
%}

else
    
y=sum(y,2);
y=y/sum(y);
plot(x2{i},y(1:end-1),sym{j},'LineWidth',1.5);


end
hold on;

  end
  
if i==2
        set(gca,'Xtick',0:2:10);
 end
    
    if i==3
        set(gca,'Xtick',0:0.4:2);
        set(gca,'XMinorTick','On');
    end


%plot(x2{i}(1:end),y2{i}(1:end-1,1:6),'k-','LineWidth',0.5);

set(gca,'FontSize',14);
xlabel(xlb{i});
ylabel('Probability');
%xlim([0,x{i}(end-2)]);
ylim(ylm(i,:));
xlim(xlm(i,:));
if i==1
     legend('WM','Subcortical');
end
end

