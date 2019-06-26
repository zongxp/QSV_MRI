function [v,sart,s_static]=calc_FlowEnhancement(TR,T1,TE,T2,thk,FA,pulse_profile)
%  [v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static]=Prep_PV_V4TOF2Sur_PC(TR,T1,TE,T2,thk,FA,pulse_profile)
% calculate the interpolation curves for PV_V4TOF2Sur_PC calculatioin.
% TR: TR of the TOF and PC scans
% T1: 1*2 vector; T1 of the flowing and static spins  in s
% TE: 
% T2: 1*2 vector; T2 of the flowing and static spins in s
% FA:  flip angle in degrees
% thk: thickness of the slice in cm
% pulse_profile: sinc or square


%%
TR=double(TR);
T1=double(T1);
T2=double(T2);
thk=double(thk);
FA=double(FA);


%%

v=[0:0.01:0.05,0.06:0.02:10,11:90];
%v_calc_pre=0.7692;
%v_calc_pre=60;
disp('precalculate the signal vs v curve');
sart=zeros(length(v),1);

for i=1:length(v)
   sart(i)=total_signal_new(v(i),TR,thk,FA,T1(1),pulse_profile,1)*exp(-TE/T2(1));
end

s_static=total_signal_new(0,TR,thk,FA,T1(2),pulse_profile,1)*exp(-TE/T2(2));  % was using total_signal; should be the same as total_signal

s0=sart(1);
sart=sart/s0;
s_static=s_static/s0;


function sart=total_signal_new(v,TR,thk,fa,T1,pulse_profile,sz_init)
% pos and fa gives the RF slice selection profile;
% for an ideal boxcar profile use pos=[-thk/2 thk/2], fa=[fa,fa];

if ~exist('sz_init','var')
    sz_init=1;
end

if strcmp(pulse_profile,'square')
    pos=[-thk/2 thk/2];
    fa=[fa,fa];
elseif strcmp(pulse_profile,'sinc')
    
        pos= linspace(-6,6,400);
        faNorm=[3.0e-05 3.0e-05 2.9e-05 2.5e-05 2.0e-05 1.4e-05 5.9e-06 2.3e-06 1.1e-05 1.9e-05 2.6e-05 3.2e-05 3.6e-05 3.8e-05 3.8e-05 3.5e-05 3.0e-05 2.3e-05 1.5e-05 4.8e-06 5.7e-06 1.6e-05 2.6e-05 3.5e-05 4.2e-05 4.7e-05 4.9e-05 4.8e-05 4.4e-05 3.7e-05 2.7e-05 1.6e-05 2.6e-06 1.1e-05 2.5e-05 3.7e-05 4.8e-05 5.7e-05 6.2e-05 6.3e-05 6.1e-05 5.5e-05 4.5e-05 3.2e-05 1.6e-05 1.7e-06 2.0e-05 3.8e-05 5.4e-05 6.8e-05 7.8e-05 8.4e-05 8.5e-05 8.0e-05 7.1e-05 5.6e-05 3.7e-05 1.5e-05 9.8e-06 3.5e-05 5.9e-05 8.1e-05 9.9e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.3e-05 7.1e-05 4.3e-05 1.0e-05 2.5e-05 6.1e-05 9.5e-05 1.2e-04 1.5e-04 1.6e-04 1.7e-04 1.7e-04 1.5e-04 1.3e-04 9.0e-05 4.6e-05 4.3e-06 5.8e-05 1.1e-04 1.6e-04 2.1e-04 2.4e-04 2.6e-04 2.7e-04 2.6e-04 2.3e-04 1.9e-04 1.3e-04 5.5e-05 2.6e-05 1.1e-04 2.0e-04 2.8e-04 3.4e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.3e-04 2.4e-04 1.2e-04 2.3e-05 1.8e-04 3.5e-04 5.1e-04 6.6e-04 7.9e-04 8.8e-04 9.3e-04 9.3e-04 8.7e-04 7.5e-04 5.7e-04 3.3e-04 5.0e-05 2.7e-04 6.0e-04 9.2e-04 1.2e-03 1.4e-03 1.6e-03 1.6e-03 1.5e-03 1.3e-03 8.7e-04 2.7e-04 4.9e-04 1.4e-03 2.4e-03 3.6e-03 4.7e-03 5.7e-03 6.6e-03 7.2e-03 7.5e-03 7.1e-03 6.1e-03 4.1e-03 1.2e-03 3.0e-03 8.6e-03 1.6e-02 2.5e-02 3.5e-02 4.8e-02 6.3e-02 8.0e-02 1.0e-01 1.2e-01 1.5e-01 1.7e-01 2.0e-01 2.3e-01 2.6e-01 3.0e-01 3.3e-01 3.7e-01 4.1e-01 4.4e-01 4.8e-01 5.2e-01 5.6e-01 6.0e-01 6.4e-01 6.7e-01 7.1e-01 7.4e-01 7.7e-01 8.0e-01 8.3e-01 8.6e-01 8.8e-01 9.0e-01 9.2e-01 9.4e-01 9.5e-01 9.7e-01 9.8e-01 9.8e-01 9.9e-01 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 9.9e-01 9.8e-01 9.8e-01 9.7e-01 9.5e-01 9.4e-01 9.2e-01 9.0e-01 8.8e-01 8.6e-01 8.3e-01 8.0e-01 7.7e-01 7.4e-01 7.1e-01 6.7e-01 6.4e-01 6.0e-01 5.6e-01 5.2e-01 4.8e-01 4.4e-01 4.1e-01 3.7e-01 3.3e-01 3.0e-01 2.6e-01 2.3e-01 2.0e-01 1.7e-01 1.5e-01 1.2e-01 1.0e-01 8.0e-02 6.3e-02 4.8e-02 3.5e-02 2.5e-02 1.6e-02 8.6e-03 3.0e-03 1.2e-03 4.1e-03 6.1e-03 7.1e-03 7.5e-03 7.2e-03 6.6e-03 5.7e-03 4.7e-03 3.6e-03 2.4e-03 1.4e-03 4.9e-04 2.7e-04 8.7e-04 1.3e-03 1.5e-03 1.6e-03 1.6e-03 1.4e-03 1.2e-03 9.2e-04 6.0e-04 2.7e-04 5.0e-05 3.3e-04 5.7e-04 7.5e-04 8.7e-04 9.3e-04 9.3e-04 8.8e-04 7.9e-04 6.6e-04 5.1e-04 3.5e-04 1.8e-04 2.3e-05 1.2e-04 2.4e-04 3.3e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.4e-04 2.8e-04 2.0e-04 1.1e-04 2.6e-05 5.5e-05 1.3e-04 1.9e-04 2.3e-04 2.6e-04 2.7e-04 2.6e-04 2.4e-04 2.1e-04 1.6e-04 1.1e-04 5.8e-05 4.3e-06 4.6e-05 9.0e-05 1.3e-04 1.5e-04 1.7e-04 1.7e-04 1.6e-04 1.5e-04 1.2e-04 9.5e-05 6.1e-05 2.5e-05 1.0e-05 4.3e-05 7.1e-05 9.3e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.9e-05 8.1e-05 5.9e-05 3.5e-05 9.8e-06 1.5e-05 3.7e-05 5.6e-05 7.1e-05 8.0e-05 8.5e-05 8.4e-05 7.8e-05 6.8e-05 5.4e-05 3.8e-05 2.0e-05 1.7e-06 1.6e-05 3.2e-05 4.5e-05 5.5e-05 6.1e-05 6.3e-05 6.2e-05 5.7e-05 4.8e-05 3.7e-05 2.5e-05 1.1e-05 2.6e-06 1.6e-05 2.7e-05 3.7e-05 4.4e-05 4.8e-05 4.9e-05 4.7e-05 4.2e-05 3.5e-05 2.6e-05 1.6e-05 5.7e-06 4.8e-06 1.5e-05 2.3e-05 3.0e-05 3.5e-05 3.8e-05 3.8e-05 3.6e-05 3.2e-05 2.6e-05 1.9e-05 1.1e-05 2.3e-06 5.9e-06 1.4e-05 2.0e-05 2.5e-05 2.9e-05 3.0e-05 3.0e-05];
        pos=pos(101:300);
        faNorm=faNorm(101:300);
    
    %faNorm=ri('SINC_Profile_Normalized4FA25.mat','','','angle');
    
    pos=pos*thk/2;
    fa=fa*faNorm;
   
end
%T1art=0.2;
%T1tissue=1.4;
pos=pos-min(pos);  % start from 0
thk2=max(pos);


if v==0
    
    fa=fa(1:4:end);
    
    for i=1:length(fa)
        sart_tmp(i)=ssFLASH(fa(i),TR,T1,0,1);
    end
    sart=mean(sart_tmp);
    
    sart=sart*thk2/thk;
    return;
end

TotalSeg=100;

step=thk2/TotalSeg;

if v*TR>step
    
    sart_tmp =zeros(1,TotalSeg);
    
    for j=1:TotalSeg
        
        z=(j-1)*step;
        
        nRF=  ceil(z/v/TR);% number of RF that has been experienced by this spin
        
        sz=sz_init;
        for irep=1:nRF
            
            z0=z-(nRF-irep)*v*TR;
            fa0=interp1(pos,fa,z0);
            
            if irep==nRF
                sart_tmp(j) = sz*sin(fa0/180*pi);
            end
            sz=sz*cos(fa0/180*pi);
            
            
            sz=1+(sz-1)*exp(-TR/T1);
        end
        
    end
    
else
    step=v*TR;
    nRF = ceil(thk2/step);
    
    sart_tmp =zeros(1,nRF);
    sz=sz_init;
    
    for j=1:nRF
        
        z=(j-1)*step;
              
        fa0=interp1(pos,fa,z);
        
        sart_tmp(j) = sz*sin(fa0/180*pi);
        sz=sz*cos(fa0/180*pi);
        
        sz=1+(sz-1)*exp(-TR/T1);
        
    end
    
end
% below is same as
sart=mean(sart_tmp(:));

sart=sart*thk2/thk;


