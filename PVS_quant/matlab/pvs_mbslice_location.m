function pvs_mbslice_location(pos1,pos2)

dist_factor=100*(abs(pos1-pos2)/2-1);


fprintf('Slice position = %3.1f\n',(pos1+pos2)/2);

fprintf('dist factor for sbRef = %d\n',round(abs(dist_factor)));

fprintf('Slice separation = %3.1f\n',abs(pos1-pos2)/2);


