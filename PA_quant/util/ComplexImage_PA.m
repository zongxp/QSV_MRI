function im=ComplexImage_PA(center,rad,vmean,fit_par)
% center in cm
% rad in mm
% vmean in cm/s
% fit_par contains:
% roi: voxels within roi will be included in fitting
% sz: matrix size of the image
% voxSize: cm
% voxSize_interp: cm
% VENC: in cm/s
% interp_factor: interpolation factor for convolution calculation
% va: array of blood velocities cm/s
% sa: signal intensity at va 


rad=rad*0.1; % to cm

sz = fit_par.sz;
voxSize=fit_par.voxSize;
voxSize_interp=fit_par.voxSize_interp;
venc=fit_par.VENC;
flow_pattern=fit_par.flow_pattern;

va=fit_par.va;  %cm/s
sa=fit_par.sa;

FOV=sz.*voxSize_interp;

fr2=@(r) fr(r,vmean,va,sa,rad,venc,flow_pattern);
im=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp,false,fit_par.interp_factor);


function res=fr(r,vmean,va,sa,rad,venc,flow_pattern)
if strcmp(flow_pattern,'Laminar')
    vi=v_laminarFlow(r/rad,vmean);
elseif strcmp(flow_pattern,'BluntedParabolic')
    vi=v_BluntedParabolicFlow(r/rad,vmean);
elseif  strcmp(flow_pattern,'Plug')
    vi=vmean*ones(size(r));
else
    error('Unknow velocity pattern');
end

si=interp1(va,sa,vi);
res=si.*(exp(1i*vi/venc*pi));
res(r>rad)=0;

