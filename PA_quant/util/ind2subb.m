
function sub=ind2subb(mtrx,ind_arr)

sub=zeros(length(ind_arr),length(mtrx));
for j=1:length(ind_arr)
     ind=ind_arr(j);
for i=1:length(mtrx)
   
    if i==length(mtrx)
        sub(j,1)=ind;
    else
    sub(j,end-i+1)=floor((ind-1)/prod(mtrx(1:end-i)))+1;
    end
    ind=ind-(sub(j,end-i+1)-1)*prod(mtrx(1:end-i));
end

end