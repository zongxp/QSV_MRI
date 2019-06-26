function curv=curvature(path,vox_size)


path=path.*repmat(vox_size,[size(path,1),1]);
path=path-mean(path,1);

d=pca(path);

%path=path-path.*d(:,3)';  % remove the third component

path = path*d;   % transform to the x-y plane;

fitfun2=@(x) fitfun(x,path);

pc=polyfit(path(:,1),path(:,2),2);
if pc(1)<0.000001
    curv=0;
else
  R = 1/abs((2*pc(1)));
  res=lsqnonlin(fitfun2,[R,0,R*sign(pc(1))]);
  curv=1/res(1);
end







function [res,d_fit]=fitfun(x,d)

rad=x(1);
center=x(2:3);
center(3)=0;
d=d-center;
phi=atan2(d(:,2),d(:,1));

d_fit=[cos(phi),sin(phi),0*phi]*rad;
res=d-d_fit;
d_fit=d_fit+center;













