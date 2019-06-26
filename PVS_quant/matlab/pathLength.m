function [len,len_all]=pathLength(path,ind,voxsize,matsize)


   pos=ind2subb(matsize,ind);

   len=0;
   len_all=zeros(1,length(path));

   for i=1:length(path)-1
    
       len=len+sos((pos(path(i+1),:)-pos(path(i),:)).*voxsize(1:3),2);
       len_all(i+1)=len;
    
   end

