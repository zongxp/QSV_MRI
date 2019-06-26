function dmax=cluster_maxdist(c)
% caculate the max extent of the cluster in c, where c is a mask with
% only one cluster.

ind=find(c>0);

s=ind2subb(size(c),ind);

dmax=0;
for i=1:size(s,1)
    for j=i+1:size(s,1)
        
        d=sos(s(i,:)-s(j,:),2);
        
        if d>dmax
            dmax=d;
        end
    end
end
    
      dmax=dmax+1;
