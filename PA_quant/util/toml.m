function d=toml


d=mfilename('fullpath');
d=fileparts(d);
if nargout==0
  cd(d);   
end