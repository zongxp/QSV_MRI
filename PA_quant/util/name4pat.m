function res=name4pat(pat)
% name4pat(pat)

dir_str=dir(pat);
dname=fileparts(pat);
if isempty(dname)
    dname=pwd;
end
if length(dir_str)>1
   % warning('more than 1 file found');
elseif isempty(dir_str)
   res=[];
   return;
end


if length(dir_str)==1
 res=fullfile(dname,dir_str(1).name);
else
    for i=1:length(dir_str)
    res{i}=fullfile(dname,dir_str(i).name);
    end
    
end

