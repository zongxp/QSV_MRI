function Fit_MACD(params,save_name)
% Fit_MACD(params,save_name)
%% params should contain the following fields:
% roi_vessel: vessel mask
% roi_wm: white matter mask
% voxSize: acquired voxSize; (cm)
% voxSize_interp: voxSize of the images; (cm)
% thk: slice thickness: (cm)
% mag: magnitude images;
% ph: phase images; (value should be in degrees).
% FA: flip angle;
% roiRad: radius for the circular ROI encompassing all pixels to be included in model image fitting (cm).
% bgRad: background ROI out- and inner- radii for calculating tissue signal (cm)
% interp_factor: interpolation factor for PSF convolution.
% VENC: (cm/s)
% T1: T1 for blood and WM tissue, respectively. (s)
% T2: T2 for blood and WM tissue, respectively. (s)
% TR: (s)
% TE: (s)
% negPhase: whether the flow produce a negative phase
% flow_pattern: 'Laminar' (default),'BluntedParabolic' or 'Plug'
% plot_fit: plot fit results
% skip_fit: skip fitting
% rad0 and v0: radius and flow velocity in mm and cm/s; initial values for
% fitting; also used for calculating the fitting image when skip_fit is
% true;
% lambda: blood-WM partition coefficient

T1=params.T1;  %blood, wm.

if ~isfield(params,'T2')
    T2=[28.5,23.5]*0.001;
else
    T2=params.T2; %  %putaman T2* is 17.6 ms
end

neg_phase=params.negPhase; %false;

thk=params.thk; % cm
TR = params.TR; % in s
TE = params.TE; %

FA=params.FA;
if ~isfield(params,'flow_pattern')
    params.flow_pattern='Laminar';
end

[va,sa,s_static]=prep_interp(TR, T1, T2, TE, thk, FA);

%%

roi_vessel=params.roi_vessel;
if max(roi_vessel(:))==1 %not clusterized
    roi_vessel=clusterize2_2d(roi_vessel);
end


nroi=max(roi_vessel(:));
v=zeros(1,nroi);
rad = zeros(1,nroi);

center=zeros(2,nroi);
mbg=zeros(2,nroi);

cd_fit=[];
data=[];

if neg_phase
    params.ph=-params.ph;
end
if nroi==0
    return;
end
for i=1:nroi
    %%
    
    tic;
    if sum(roi_vessel(:)==i)==0
        continue;
    end
    
    
    [data(:,:,i),roi,mbg(:,i),data_debug(i)]=calc_cd(params,roi_vessel,i,s_static);
    
    fit_par=params;
    fit_par.va=va;
    fit_par.sa=sa;
    fit_par.data=data(:,:,i);
    fit_par.roi=roi;
    fit_par.sz=size(roi);
    fitfunc2 = @(x) MBAC_fitfunc(x,fit_par);
    options=prep_options;
    
    % x(1:2) in cm;  x(3) in mm;  x(4) in 10 cm/s
    if isfield(params,'skip_fit') && params.skip_fit
        res=[0,0,params.rad0,params.v0/10];
    else
        [res,tmp,tmp2,exitflag(i),output(i)]=lsqnonlin(fitfunc2,[0,0,params.rad0,params.v0/10],[-0.05,-0.05,0.01,0],[0.05,0.05,0.25,0.5],options);
    end
    
    center(:,i) = res(1:2)*10; % in mm
    v(i) = res(4)*10;  % in cm/s
    rad(i)=res(3); %% in mm
    
    [resid,cd_fit(:,:,i)]=fitfunc2(res);
    
    if isfield(params,'plot_fit') && params.plot_fit  %check fit
        plot_fit(data(:,:,i),cd_fit(:,:,i),roi);
    end
    fprintf('vessel %d/%d: v = %s cm/s; d = %s mm\n',i,nroi,num2str(v(i)),num2str(rad(i)*2));
    fprintf('Time remaining %4.3f s\n',toc*(nroi-i));
end

roi_vessel=params.roi_vessel;
params=rmfield(params,{'mag','ph','roi_vessel'});

save(save_name,'v','rad','cd_fit','data','center','params','exitflag','output','roi_vessel','data_debug');

function options=prep_options
if datenum(version('-date'))>=736580
    options = optimoptions('lsqnonlin','MaxIter',Inf,'StepTolerance',1e-6,...
        'MaxFunctionEvaluations',Inf,'FunctionTolerance',0,...
        'OptimalityTolerance',0,'Display','iter',...
        'FiniteDifferenceStepSize',0.01);
else
    options = optimoptions('lsqnonlin','MaxIter',Inf,'TolX',1e-6,...
        'MaxFunEvals',Inf,'TolFun',0,'TolPCG',0,'Display','iter','FinDiffRelStep',0.01);
end

function plot_fit(data,cd_fit,roi)

data(~roi)=0;

cd_fit(~roi)=0;
figure; compare_cd(data,cd_fit);
resid=sos(data(:)-cd_fit(:),1).^2;
disp(resid);

function [cd,roi3,mbg,data_debug]=calc_cd(p,roi_vessel,i,s_static)

roi2=roi_vessel==i;
p.bgRad=sort(p.bgRad,'descend');
bgRad_i=round(p.bgRad./p.voxSize_interp);

bg_out=mask_circle(size(roi2),bgRad_i(1),roiCOM(roi2),1);
bg_in=mask_circle(size(roi2),bgRad_i(2),roiCOM(roi2),1);

bg=bg_out&~bg_in;

bgroi=bg&p.roi_wm(:,:)>0&roi_vessel==0;

if sum(bgroi(:))<118 % if partly outside wm, do not pose restriction.
    bgroi=bg;%
end

pos=roiCOM(roi2);
roiRad_i=round(mean(p.roiRad./p.voxSize_interp));  % 1 and 2 should be equal

cropr=pos(1)-roiRad_i:pos(1)+roiRad_i;
cropc=pos(2)-roiRad_i:pos(2)+roiRad_i;

cropr_l=pos(1)-bgRad_i(1)-5:pos(1)+bgRad_i(1)+5;
cropc_l=pos(2)-bgRad_i(1)-5:pos(2)+bgRad_i(1)+5;

cropr_l(cropr_l<1)=1;
cropr_l(cropr_l>size(bg,1))=size(bg,1);

cropc_l(cropc_l<1)=1;
cropc_l(cropc_l>size(bg,2))=size(bg,2);


if 1
    if size(p.ph,4)>1 % phase images
        ph=tof_bg_rm(p.ph(:,:,1,2*i-1:2*i),bgroi,[],[],false);
        mag=tof_bg_rm(p.mag(:,:,1,2*i-1:2*i),bgroi);
    else % phase difference images
        ph=tof_bg_rm(p.ph(:,:),bgroi,[],[],false);
        mag=tof_bg_rm(p.mag(:,:,1,:),bgroi);
    end
else
    mag=p.mag;
    ph=p.ph;
end

mbg=mean_roi(mag,bgroi>0);

snorm=mean(mbg)/s_static*p.lambda;

mag1=mag(:,:,:,1);
mag2=mag(:,:,:,2);

if size(ph,4)==1
    cd=mag2.*exp(1i*ph/180*pi) - mag1;
else
    cd=mag2.*exp(-1i*ph(:,:,:,2)/180*pi) - mag1.*exp(-1i*ph(:,:,:,1)/180*pi);
end

cd=cd/snorm;

%         z(1)=mean(cd(roi2));
%         z(2)=mean(mag1(roi2))/snorm;
%

%% debug purpose
data_debug.ph_crop=ph(cropr_l,cropc_l,:);
data_debug.mag1_crop=mag(cropr_l,cropc_l,:,1);
data_debug.mag2_crop=mag(cropr_l,cropc_l,:,2);
data_debug.cd_crop=cd(cropr_l,cropc_l);
data_debug.bgroi_crop=bgroi(cropr_l,cropc_l);
roi3=mask_circle(size(cd),roiRad_i,pos,1);
data_debug.roi_crop=roi3(cropr_l,cropc_l);


roi3=roi3(cropr,cropc);
cd=cd(cropr,cropc);






