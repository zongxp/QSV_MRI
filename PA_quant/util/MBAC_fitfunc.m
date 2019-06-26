function [res,cd_fit]=MBAC_fitfunc(x,fit_par)
% x(1:2) in cm
% x(3) in mm
% x(4) in 10 cm/s
% fit_par contains:
    % roi: voxels within roi will be included in fitting
    % data: CD image
    % voxSize: cm
    % voxSize_interp: cm
    % VENC: in cm/s
    % interp_factor: interpolation factor for convolution calculation
    % flow_pattern: Plug, Lamina, or BluntedParobolic

cd=fit_par.data;
roi= fit_par.roi;

im1=ComplexImage_PA(x(1:2),x(3),x(4)*10,fit_par);

fit_par.VENC=Inf;
im0=ComplexImage_PA(x(1:2),x(3),x(4)*10,fit_par);

cd_fit=im1-im0;

resid=cd(roi>0)-cd_fit(roi>0);
res=[real(resid(:));imag(resid(:))];

