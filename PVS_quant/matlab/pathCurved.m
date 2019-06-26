function pathCurved(fname)
% file should contain a ind file (1*npath)
% path file (1*npath)
% c: non zero with values equal to the npath index.

a=load(fname);
ind=a.ind;
sz=size(a.c);
nfound=0;
nchecked=0;
path=a.path;
fixedpath=zeros(sz);
probclust=fixedpath;

while nchecked<length(path)
   
    if mod(nchecked+1,50)==0
    %    disp(nchecked);
    end
    
    i=nchecked+1;
    
    xpath=setdiff(1:length(ind{i}),path{i});
     
       ind{i}=ind{i}(:);
    ind{i}=[ind{i}(path{i});ind{i}(xpath)];
 
    path{i}=1:length(path{i});
    
   pospath=ind2subb(sz,ind{i}(path{i}));
    
   breakpath=false;
   for j=1:size(pospath,1)-3
       
       d1=pospath(j,:)-pospath(j+1,:);
       d2=pospath(j+2,:)-pospath(j+3,:);

   
       cs=sum(d1.*d2)/sos(d1)*sos(d2);
       
       if cs<0
           probclust(ind{i})=1;
           [indnew,pathnew]=splitPath(ind{i},path{i},[j+1,j+2],sz);
           
           path(i)=[];
           ind(i)=[];
           path(end+1:end+2)=pathnew;
           ind(end+1:end+2)=indnew;
           fixedpath(indnew{end})=1;
           fixedpath(indnew{end-1})=1;
           
           nfound=nfound+1;
           breakpath=true;
           break;
           
       end
       
   end
   
   if breakpath
       continue;
   end
   
   
   
   for j=1:size(pospath,1)-4
       
       d1=pospath(j,:)-pospath(j+1,:);
       d2=pospath(j+3,:)-pospath(j+4,:);
   
       cs=sum(d1.*d2)/sos(d1)*sos(d2);
     
       if cs<0
            probclust(ind{i})=1;
           [indnew,pathnew]=splitPath(ind{i},path{i},j+2,sz);
           path(i)=[];
           ind(i)=[];
           path(end+1:end+2)=pathnew;
           ind(end+1:end+2)=indnew;
           
           fixedpath(indnew{end})=1;
           fixedpath(indnew{end-1})=1;
           
           nfound=nfound+1;
           breakpath=true;
           break;
       end
   end
   
   
   if breakpath
       continue;
   end
   
   nchecked=nchecked+1;
   
end

if nfound>0
fprintf('%d paths found with greater than 90 degree curvature\n',nfound);
end

fname=strtok(fname,'.');

c=zeros(sz);
pathroi=zeros(sz);

for i=1:length(path)

    c(ind{i})=i;
    
    pathroi(ind{i}(path{i}))=i;

end


save([fname,'_Curv'],'c','ind','pathroi','path','fixedpath','probclust');


function [indnew,pathnew]=splitPath(ind,p,n,sz)


pathnew=cell(1,2);
pathnew{1}=p(1:n(1)-1);
pathnew{2}=p(n(end)+1:end);


indnew{1}=ind(pathnew{1});
indnew{2}=ind(pathnew{2});
pathnew{1}=1:length(pathnew{1});
pathnew{2}=1:length(pathnew{2});
indpath1=indnew{1};
indpath2=indnew{2};

ind(p(n))=[];

ind3=setdiff(ind,[indnew{1};indnew{2}]);

for i=1:size(ind3,1)
    
    
    if minDist(ind3(i),indpath1,sz)<minDist(ind3(i),indpath2,sz)
      indnew{1}(end+1)=ind3(i);
    else
      indnew{2}(end+1)=ind3(i);
    end

       
    
end

    


    function res=minDist(ind0,pathind,sz)
        
          pos1=ind2subb(sz,ind0);
          pos2=ind2subb(sz,pathind);
            
          dist=sos(repmat(pos1,[size(pos2,1),1])-pos2);
            
          res=min(dist);
                
                
        


