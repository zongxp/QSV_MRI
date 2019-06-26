function vloc=v_BluntedParabolicFlow(r,vmean)

% equation (9) in Koutsiaris Clinical Hemorheology and Microcirculation 43 (2009) 321–334 

vmax=vmean*1.4897;

vloc=vmax*(1-0.58*r.^2).*(1-r.^22);

vloc(r>1)=0;

% %% 1.4897 calculated as below; 
% f=@(r) (1-0.58*r.^2).*(1-r.^22).*2.*r;
% a=integral(f,0,1);
% disp(1/a);
