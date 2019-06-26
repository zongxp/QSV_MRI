function [res,res0]=circleImage_CartesianKSpace(fr,center,FOV,voxSize,voxSize_interp,angular_dependent,interp)
% res=circleImage_CartesianKSpace(fr,center,FOV,voxSize,voxSize_interp)
% fr: the function of f(radius) or f(radius,phi) where phi is measured from
% x axis which is the second dimension.
% center: the center of the circle; mm
% FOV: field of view; mm
% voxSize: acquired voxel size; in mm
% voxSize_interp: voxel size for output res; in mm
% interp: interpolation factor for calculating the PSF; relative to
% voxSize_interp
% 
if ~exist('angular_dependent','var')
    angular_dependent = false;
end

if ~exist('interp','var')
    interp=10;
end

if length(FOV)==1
    FOV=FOV*[1,1];
end

if length(voxSize)==1
    voxSize=voxSize*[1,1];
end

if length(voxSize_interp)==1
    voxSize_interp=voxSize_interp*[1,1];
end

nacq=FOV./voxSize;
dk=1./FOV*2*pi;

nvox=round(nacq.*voxSize./voxSize_interp);

%kmax=dk.*nacq/2;
kmax=pi./voxSize;


res=zeros(nvox);

%dxy_interp=dxy/interp;

% xk=linspace(-lim(1),lim(1),round(lim(1)/dxy(1)*5));
% yk=linspace(-lim(2),lim(2),round(lim(2)/dxy(2)*5));

x0=linspace(-FOV(2)/2,FOV(2)/2,nvox(2)*interp+1);
y0=linspace(-FOV(1)/2,FOV(1)/2,nvox(1)*interp+1);
x0=x0(1:end-1);
y0=y0(1:end-1);

x=repmat(x0,[length(y0),1]);
y=repmat(y0',[1,length(x0)]);
% 
% xk=x0(x0>=-lim(1)&x0<lim(1));
% yk=y0(y0>=-lim(2)&y0<lim(2));

dxy=FOV./nvox/interp;

lim=3*pi./kmax;%3*voxsize

xk=linspace(-lim(2),lim(2),2*lim(2)/dxy(2)+1);
xk=xk(1:end-1);
yk=linspace(-lim(1),lim(1),2*lim(1)/dxy(1)+1);
yk=yk(1:end-1);

nxk=length(xk);
nyk=length(yk);

xkm=repmat(xk,[length(yk),1]);% x is row
ykm=repmat(yk',[1,length(xk)]);

kernelx=kernelSinc(kmax(2),xkm);
kernely=kernelSinc(kmax(1),ykm);

sumx=sum(kernelx(:))/size(kernelx,1);
sumy=sum(kernely(:))/size(kernely,2);

res0=fr2(fr,x-center(2),y-center(1),angular_dependent);
nx=size(res0,2);
ny=size(res0,1);


for i=1:nvox(1)
    for j=1:nvox(2)
        

       
     p(1)=(i-half(nvox(1)))*interp;
     p(2)=(j-half(nvox(2)))*interp;
 
     ind1=(1:ny)+p(1)+half(nyk)-half(ny);
     ind2=(1:nx)+p(2)+half(nxk)-half(nx);
     
     sel1=ind1>1 & ind1<nyk;
     sel2=ind2>1 & ind2<nxk;
     
     frtmp=res0(sel1,sel2);
     
     tmp=kernelx(ind1(sel1),ind2(sel2)).*kernely(ind1(sel1),ind2(sel2)).*frtmp;
     res(i,j)=sum(tmp(:))/sumx/sumy;
     
      
    end
    
end

function res=half(x)

   if mod(x,2) == 0
       res=x/2+1;
   else
      res = (x+1)/2; 
   end
       

function res=fr2(fr,x,y,angular_dependent)

if angular_dependent
   phi= atan2(y,x)*180/pi;
   r=sqrt(x.^2+y.^2);
   res=fr(r,phi);
    
else
    
r=sqrt(x.^2+y.^2);
res=fr(r);
end
% 
% function res=fr2(rad,x,y)
% 
% res=ones(size(x));
% r=sqrt(x.^2+y.^2);
% res(r>rad)=0;




