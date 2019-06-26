function oldstr = strrm(oldstr,substr)
%
% newstr = strrm(oldstr,substr)
%remove substring substr from the string oldstr.

while 1
 i=strfind(oldstr,substr);
 if isempty(i)
     break;
 end
 oldstr(i:i+length(substr)-1)=[];
end





