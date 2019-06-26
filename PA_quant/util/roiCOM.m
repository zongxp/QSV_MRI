function res=roiCOM(m)
 
ind=ind2subb(size(m),find(m>0));

if ~isempty(ind)
res=round(mean(ind,1));  % center of mass

else
    res=0;
end

