function pathCurved_nofix(fname)
% file should contain a ind variable (1*npath)
% i_ind_path variable (1*npath)
% c: non zero with values equal to the npath index.

a=load(fname);
ind=a.ind;
sz=size(a.c);
nfound=0;
nchecked=0;
i_ind_path=a.i_ind_path;

probclust=zeros(sz);
niter=0;
while nchecked<length(i_ind_path)
   
    niter=niter+1;
    if mod(nchecked+1,50)==0
    %    disp(nchecked);
    end
    
    i=nchecked+1;
    
    xpath=setdiff(1:length(ind{i}),i_ind_path{i});
     
       ind{i}=ind{i}(:);
    ind{i}=[ind{i}(i_ind_path{i});ind{i}(xpath)];
 
    i_ind_path{i}=1:length(i_ind_path{i});
    
   pospath=ind2subb(sz,ind{i}(i_ind_path{i}));
    
   breakpath=false;
   for j=1:size(pospath,1)-3
       
       d1=pospath(j,:)-pospath(j+1,:);
       d2=pospath(j+2,:)-pospath(j+3,:);

   
       cs=sum(d1.*d2)/sos(d1)/sos(d2);
       
       if cs<0
           d1=pospath(j,:)-pospath(j+1,:);
           d3=pospath(j+1,:)-pospath(j+2,:);      
           cs2=sum(d1.*d3)/sos(d1)/sos(d3);
           if cs2<0
               
               probclust(ind{i})=1;
               %  [indnew,pathnew]=splitPath(ind{i},path{i},[j+1,j+2],sz);
               
               i_ind_path(i)=[];
               ind(i)=[];
               %  path(end+1:end+2)=pathnew;
               %  ind(end+1:end+2)=indnew;
               %  fixedpath(indnew{end})=1;
               %  fixedpath(indnew{end-1})=1;
               disp(pospath(j,:));
               disp(pospath(j+1,:));
               disp(pospath(j+2,:));
               
               nfound=nfound+1;
               fprintf('Path %d angle > 90o; removed; angle = %f\n',i,acos(cs)*180/pi);
               breakpath=true;
               break;
           end
       end
       
   end
   
   if breakpath
       continue;
   end
   
   
   %{
   for j=1:size(pospath,1)-4
       
       d1=pospath(j,:)-pospath(j+1,:);
       d2=pospath(j+3,:)-pospath(j+4,:);
   
       cs=sum(d1.*d2)/sos(d1)/sos(d2);
     
       if cs<0
            probclust(ind{i})=1;
         %  [indnew,pathnew]=splitPath(ind{i},path{i},j+2,sz);
           path(i)=[];
           ind(i)=[];
          % path(end+1:end+2)=pathnew;
          % ind(end+1:end+2)=indnew;
           
          % fixedpath(indnew{end})=1;
          % fixedpath(indnew{end-1})=1;
           
           nfound=nfound+1;
           breakpath=true;
           break;
       end
   end
   
   %}
   if breakpath
       continue;
   end
   
   nchecked=nchecked+1;
   
end

if nfound>0
fprintf('%d paths found with >90 degree angle\n',nfound);

else
    disp('No paths with > 90 degree angle found.  Good!');
end

fname=strtok(fname,'.');

c=zeros(sz);
pathroi=zeros(sz);

for i=1:length(i_ind_path)

    c(ind{i})=i;
    
    pathroi(ind{i}(i_ind_path{i}))=i;

end

voxsize=a.voxsize;
save([fname,'_Curv'],'c','ind','pathroi','i_ind_path','probclust','voxsize');


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
                
                
        


