function res=param2name(varargin)

res=num2str(varargin{1});

for i=2:nargin
    
    res=[res,'_',num2str(varargin{i})];
    
end

res=strrep(res,' ','_');
res=rm_rep_ch(res,'_');
res=[res,'.mat'];

