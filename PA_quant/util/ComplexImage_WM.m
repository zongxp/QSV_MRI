function im=ComplexImage_WM(center,rad,fit_par)
% center in cm
% rad in mm of PA within WM
% fit_par contains:
% sz: matrix size of the image
% voxSize: cm
% voxSize_interp: cm
% interp_factor: interpolation factor for convolution calculation
% st: tissue signal;


rad=rad*0.1; % to cm

voxSize=fit_par.voxSize;
voxSize_interp=fit_par.voxSize_interp;

FOV=fit_par.sz.*voxSize_interp;

ftrue1=@(x,y) signal_static(x-center(1),y-center(2),rad,fit_par.st);          
im=Image_CartesianKSpace_ConvWtFFT(ftrue1,FOV,voxSize,voxSize_interp,fit_par.interp_factor);


function res=signal_static(x,y,rad,s_static)  % do it with interp
      
    r=sqrt(x.^2+y.^2);
   
    res=0*x;
    res(r>=rad)=s_static;
    

    
          