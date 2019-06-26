function [res,ph_fit]=MAP_fitfunc(x,fit_par)
% x(1:2) in cm
% x(3) in mm
% x(4) in 10 cm/s
% fit_par contains:
    % roi: voxels within roi will be included in fitting
    % data: ph image
    % voxSize: cm
    % voxSize_interp: cm
    % VENC: in cm/s
    % interp_factor: interpolation factor for convolution calculation
    % flow_pattern: Plug, Lamina, or BluntedParobolic

ph=fit_par.data;
roi= fit_par.roi;
lambda=fit_par.lambda;

im_wm=ComplexImage_WM(x(1:2),x(3),fit_par);
im1_pa=ComplexImage_PA(x(1:2),x(3),x(4)*10,fit_par);
fit_par.VENC=Inf;
im0_pa=ComplexImage_PA(x(1:2),x(3),x(4)*10,fit_par);


im1=im1_pa*lambda+im_wm;
im0=im0_pa*lambda+im_wm;

ph_fit=angle(conj(im0).*im1)*180/pi;

resid=ph(roi>0)-ph_fit(roi>0);
res=[real(resid(:));imag(resid(:))];

