function res=tortuosity(path,vox_size)

path=path.*repmat(vox_size,[size(path,1),1]);

dp=path(1:end-1,:)-path(2:end,:);
dist=sum(sos(dp,2));

dp2=path(1,:)-path(end,:);
dist2=sos(dp2,2);

res=dist/dist2;



