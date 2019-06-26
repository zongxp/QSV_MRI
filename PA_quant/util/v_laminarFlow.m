function vloc=v_laminarFlow(r,vmean)

% determine the local velocity at radius r (normalized by tube radius)

vloc=vmean*2*(1-r.^2);  %
vloc(r>1)=0;
