function s=ssFLASH(fa,TR,T1,TE,T2s)
%s=ssFLASH(fa,TR[,T1,TE,T2s])
% fa in degrees.
% s is already the transverse magnetization after flip

if ~exist('T1','var')
    T1=2;
end

if ~exist('TE','var')
    TE=0;
end

if  ~exist('T2s','var')
   T2s = 0.025;
end
fa=fa*pi/180;
e1=exp(-TR./T1);
e2=exp(-TE./T2s);
s=sin(fa).*(1-e1)./(1-cos(fa)*e1)*e2;  %% already the transverse magnetization after flip

