function sc=mask_subcortical(roi90_fname,roi90_dir)
   cur_dir=pwd;
if isa(roi90_fname,'char')
 
    cd(roi90_dir);
    prefix=strtok(roi90_fname,'.');
    a=ri(roi90_fname);
else
    a=roi90_fname;
end

if max(a(:))==90
m=a>=71&a<=78;
else
   ml= a>=9 & a<=13; % left thalamus, thalamus_proper, caudate, palladium, putaman
   mr= a>=48 & a<=52;% right thalamus, thalamus_proper, caudate, palladium, putaman
end

[sc_l]=mask_sc(ml);
[sc_r]=mask_sc(mr);

sc_l(sc_r>0 & sc_l>0)=0;
sc_r(sc_r>0 & sc_l>0)=0;

sc=sc_l+sc_r*2;

if nargout==0
 save ([prefix,'_subcortical'], 'sc');
end


cd(cur_dir);


function [sc,ic]=mask_sc(ml)


ml2=bwmorph3d(ml,'dilate',15);
ml3=bwmorph3d(ml2,'shrink',10);
sc=ml3|ml;
ic=ml3&~ml;
