function m=ri_d1(fname,varargin)
%ri(fname,ns,do_sort,field_name)
 % ns: number of slices in the image; when empty; will be set to the
 % maximum possible.
% var_name: used for reading .mat file; which variable to read, default:roi 

    try
        
        if length(fname)<4 || ~strcmp(fname(end-3:end),'.mat')
            fname2=[fname,'.mat'];
        else
            fname2=fname;
        end
        
        m=ri(fname2,'','','d');
    catch
        m=ri(fname,varargin{:});
    end
    